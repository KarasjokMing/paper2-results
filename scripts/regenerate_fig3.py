import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from config import CHANNELS, PARAM_NAMES
from physics import compute_numpy


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT_DIR = os.path.join(ROOT, "figures")

NIGHTFIRE = LinearSegmentedColormap.from_list(
    "nightfire",
    [
        (0.00, "#02020a"),
        (0.20, "#160016"),
        (0.42, "#3c001e"),
        (0.66, "#8c0718"),
        (0.84, "#e0401a"),
        (1.00, "#fff2a6"),
    ],
)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "figure.dpi": 180,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.7,
})


def base_vector(channel):
    if channel == "T2":
        vals = {
            "Gamma": 1.0,
            "Delta": 0.008604887811638484,
            "theta_a1": 0.5020525510717644,
            "theta_a2": 1.0,
            "theta_b2": 0.0013695404362945585,
            "theta_b3": 1.0,
            "phi1": 1.502052163525476,
            "phi2": 0.9986304834268034,
            "phi3": 0.0,
        }
    else:
        vals = {
            "Gamma": 1.0,
            "Delta": 0.008088580004357304,
            "theta_a1": 0.501931278715504,
            "theta_a2": 1.0,
            "theta_b2": 1.9987126395436206,
            "theta_b3": 1.0,
            "phi1": 1.5019300301245533,
            "phi2": 0.9987126574102942,
            "phi3": 0.0,
        }
    return np.array([vals[name] for name in PARAM_NAMES], dtype=float)


def c2_vector(channel, q, eta):
    theta_b2 = q + eta if channel == "T2" else 2.0 - q + eta
    vals = {
        "Gamma": 1.0,
        "Delta": 2.0 * np.sin(np.pi * q),
        "theta_a1": 0.5 + 1.5 * q,
        "theta_a2": 1.0,
        "theta_b2": theta_b2,
        "theta_b3": 1.0,
        "phi1": 1.5 + 1.5 * q,
        "phi2": 1.0 - q,
        "phi3": 0.0,
    }
    return np.array([vals[name] for name in PARAM_NAMES], dtype=float)


def loose_vector(channel, q, eta):
    vec = base_vector(channel)
    vec[PARAM_NAMES.index("phi2")] = 1.0 - q
    if channel == "T2":
        vec[PARAM_NAMES.index("theta_b2")] = q + eta
    else:
        vec[PARAM_NAMES.index("theta_b2")] = 2.0 - q + eta
    return vec


def branch_center(channel):
    base = base_vector(channel)
    q0 = 1.0 - base[PARAM_NAMES.index("phi2")]
    if channel == "T2":
        eta0 = base[PARAM_NAMES.index("theta_b2")] - q0
    else:
        eta0 = base[PARAM_NAMES.index("theta_b2")] - (2.0 - q0)
    return q0, eta0


def heatmap(channel, constrained, q_range=(0.0005, 0.06), eta_range=(-0.02, 0.02), n=320):
    ci = CHANNELS.index(channel)
    qs = np.linspace(q_range[0], q_range[1], n)
    etas = np.linspace(eta_range[0], eta_range[1], n)
    q0, eta0 = branch_center(channel)
    if q_range[0] <= q0 <= q_range[1]:
        qs[np.argmin(np.abs(qs - q0))] = q0
    if eta_range[0] <= eta0 <= eta_range[1]:
        etas[np.argmin(np.abs(etas - eta0))] = eta0
    grid = np.empty((n, n), dtype=float)
    for iy, eta in enumerate(etas):
        for ix, q in enumerate(qs):
            vec = c2_vector(channel, q, eta) if constrained else loose_vector(channel, q, eta)
            p = compute_numpy(vec[0], vec[1], *vec[2:])
            grid[iy, ix] = p[ci] if np.all(np.isfinite(p)) else np.nan
    return qs, etas, grid


def contour_levels(grid):
    levels = [0.9, 0.99, 0.999]
    lo, hi = np.nanmin(grid), np.nanmax(grid)
    return [level for level in levels if lo <= level <= hi]


def draw_contours(ax, qs, etas, grid):
    visible = contour_levels(grid)
    if visible:
        color_map = {0.9: "#ffffff", 0.99: "#ffd166", 0.999: "#7dd3fc"}
        width_map = {0.9: 0.45, 0.99: 0.55, 0.999: 0.7}
        ax.contour(
            qs,
            etas,
            grid,
            levels=visible,
            colors=[color_map[level] for level in visible],
            linewidths=[width_map[level] for level in visible],
        )


def area_text(area):
    if area == 0:
        return "0"
    if area < 0.005:
        return f"{area:.1e}"
    return f"{area:.2f}"


def add_unconstrained_inset(ax, channel):
    q0, eta0 = branch_center(channel)
    q_half = 0.0007
    eta_half = 0.0018
    qlo = max(0.0005, q0 - q_half)
    qhi = min(0.06, q0 + q_half)
    elo = max(-0.02, eta0 - eta_half)
    ehi = min(0.02, eta0 + eta_half)
    zqs, zetas, zgrid = heatmap(channel, False, q_range=(qlo, qhi), eta_range=(elo, ehi), n=220)

    ax.add_patch(
        Rectangle((qlo, elo), qhi - qlo, ehi - elo, fill=False, edgecolor="white", linewidth=0.65)
    )
    inset = ax.inset_axes([0.58, 0.52, 0.36, 0.39])
    inset.pcolormesh(zqs, zetas, zgrid, shading="auto", cmap=NIGHTFIRE, vmin=0, vmax=1, rasterized=True)
    draw_contours(inset, zqs, zetas, zgrid)
    inset.set_xlim(qlo, qhi)
    inset.set_ylim(elo, ehi)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.text(0.06, 0.9, "zoom", transform=inset.transAxes, ha="left", va="top", fontsize=6, color="white")
    for spine in inset.spines.values():
        spine.set_edgecolor("white")
        spine.set_linewidth(0.65)


fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), constrained_layout=True)
panels = [
    ("T2", False, r"$T_2$: only $(q,\eta)$ scanned"),
    ("T2", True, r"$T_2$: $\mathcal{C}_2(q)$ enforced"),
    ("R2", False, r"$R_2$: only $(q,\eta)$ scanned"),
    ("R2", True, r"$R_2$: $\mathcal{C}_2(q)$ enforced"),
]

last_im = None
for ax, (channel, constrained, title) in zip(axes.ravel(), panels):
    qs, etas, grid = heatmap(channel, constrained)
    last_im = ax.pcolormesh(qs, etas, grid, shading="auto", cmap=NIGHTFIRE, vmin=0, vmax=1, rasterized=True)
    draw_contours(ax, qs, etas, grid)
    if not constrained:
        add_unconstrained_inset(ax, channel)
    area99 = np.nanmean(grid > 0.99)
    area999 = np.nanmean(grid > 0.999)
    ax.text(0.98, 0.05, rf"$A_{{.99}}={area_text(area99)}$" + "\n" + rf"$A_{{.999}}={area_text(area999)}$",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
            color="white", bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.55, edgecolor="none"))
    ax.set_title(title)
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$\eta/\pi$")

cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), shrink=0.92, pad=0.015)
cbar.set_label("target probability")

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT_DIR, f"fig3_correlated_branch.{ext}"))
plt.close(fig)
