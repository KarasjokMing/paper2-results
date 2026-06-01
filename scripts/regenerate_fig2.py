import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from config import CHANNELS, PARAM_NAMES
from physics import compute_numpy


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PARAM_PATH = os.path.join(ROOT, "data", "optimized_params.json")
OUT_DIR = os.path.join(ROOT, "figures")

with open(PARAM_PATH, "r", encoding="utf-8") as f:
    OPT = json.load(f)

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

LABEL = {
    "phi1": r"$\varphi_1/\pi$",
    "phi2": r"$\varphi_2/\pi$",
    "theta_a1": r"$\theta_{a1}/\pi$",
    "theta_b2": r"$\theta_{b2}/\pi$",
}

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


def vector_for(channel):
    return np.array([OPT[channel][name] for name in PARAM_NAMES], dtype=float)


def phase_map(channel, xname, yname, xr, yr, n=260):
    ci = CHANNELS.index(channel)
    xi = PARAM_NAMES.index(xname)
    yi = PARAM_NAMES.index(yname)
    base = vector_for(channel)
    xs = np.linspace(xr[0], xr[1], n)
    ys = np.linspace(yr[0], yr[1], n)
    xo = base[xi]
    yo = base[yi]
    if xr[0] <= xo <= xr[1]:
        xs[np.argmin(np.abs(xs - xo))] = xo
    if yr[0] <= yo <= yr[1]:
        ys[np.argmin(np.abs(ys - yo))] = yo
    grid = np.empty((n, n), dtype=float)
    for iy, y in enumerate(ys):
        vec = base.copy()
        vec[yi] = y
        for ix, x in enumerate(xs):
            vec[xi] = x
            p = compute_numpy(vec[0], vec[1], *vec[2:])
            grid[iy, ix] = p[ci] if np.all(np.isfinite(p)) else np.nan
    return xs, ys, grid, xo, yo


PLOTS_ALL6 = [
    ("T1", "phi1", "theta_a1", (0, 2), (0, 2), r"$T_1$: full phase plane"),
    ("R1", "phi2", "theta_a1", (-1, 1), (0, 2), r"$R_1$: periodic seam centered"),
    ("T2", "phi1", "theta_a1", (0.5, 2.0), (0.0, 1.25), r"$T_2$: extended branch view"),
    ("R2", "phi1", "theta_a1", (0.5, 2.0), (0.0, 1.25), r"$R_2$: extended branch view"),
    ("T3", "phi2", "theta_b2", (0.96, 1.04), (0.96, 1.04), r"$T_3$: near-resonant island"),
    ("R3", "phi2", "theta_b2", (0.96, 1.04), (0.96, 1.04), r"$R_3$: near-resonant island"),
]

PLOTS_MAIN = [PLOTS_ALL6[i] for i in (0, 1, 2, 4)]


def draw_maps(plots, shape, figsize, basename):
    fig, axes = plt.subplots(*shape, figsize=figsize, constrained_layout=True)
    axes_flat = np.atleast_1d(axes).ravel()
    last_im = None
    for ax, (ch, xname, yname, xr, yr, title) in zip(axes_flat, plots):
        n = 360 if ch in ("T3", "R3") else 280
        xs, ys, grid, xo, yo = phase_map(ch, xname, yname, xr, yr, n=n)
        im = ax.pcolormesh(xs, ys, grid, shading="auto", cmap=NIGHTFIRE, vmin=0.0, vmax=1.0, rasterized=True)
        last_im = im
        levels = [0.9, 0.99, 0.999]
        finite = np.isfinite(grid)
        if finite.any():
            lo, hi = np.nanmin(grid), np.nanmax(grid)
            visible_levels = [level for level in levels if lo <= level <= hi]
            if visible_levels:
                colors = ["#ffffff", "#ffd166", "#7dd3fc"][:len(visible_levels)]
                ax.contour(xs, ys, grid, levels=visible_levels, colors=colors, linewidths=[0.45, 0.55, 0.7][:len(visible_levels)])
        ax.set_title(title)
        ax.set_xlabel(LABEL[xname])
        ax.set_ylabel(LABEL[yname])
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        ax.set_xticks([xr[0], sum(xr) / 2, xr[1]])
        ax.set_yticks([yr[0], sum(yr) / 2, yr[1]])
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))

    for ax in axes_flat[len(plots):]:
        ax.axis("off")

    cbar = fig.colorbar(last_im, ax=axes_flat[:len(plots)].tolist(), shrink=0.92, pad=0.015)
    cbar.set_label("target probability")
    cbar.ax.tick_params(labelsize=7)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT_DIR, f"{basename}.{ext}"))
    plt.close(fig)


draw_maps(PLOTS_MAIN, (2, 2), (6.2, 5.0), "fig2_representative_maps")
draw_maps(PLOTS_ALL6, (2, 3), (7.4, 5.0), "fig2_all6_supplementary")
