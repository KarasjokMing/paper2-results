import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT_DIR = os.path.join(ROOT, "figures")
REPORT_DIR = os.path.join(ROOT, "reports")
PARAM_PATH = os.path.join(ROOT, "data", "optimized_params.json")

sys.path.insert(0, os.path.join(ROOT, "src"))
from direct_oracle import build  # noqa: E402


PARAM_NAMES = [
    "Gamma", "Delta", "theta_a1", "theta_a2", "theta_b2",
    "theta_b3", "phi1", "phi2", "phi3",
]
CHANNEL_TO_ROW = {"R1": 0, "T1": 1, "R2": 2, "T2": 3, "R3": 4, "T3": 5}
PI = np.pi

with open(PARAM_PATH, "r", encoding="utf-8") as f:
    OPT = json.load(f)

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


def vector_for(channel):
    return np.array([OPT[channel][name] for name in PARAM_NAMES], dtype=float)


def target_probability(channel, vec):
    args = vec.copy()
    args[2:] *= PI
    A, b = build(*args, port=1)
    out = np.linalg.solve(A, b)[6:12]
    return float(abs(out[CHANNEL_TO_ROW[channel]]) ** 2)


def c2_vector(q, eta):
    vals = {
        "Gamma": 1.0,
        "Delta": 2.0 * np.sin(np.pi * q),
        "theta_a1": 0.5 + 1.5 * q,
        "theta_a2": 1.0,
        "theta_b2": q + eta,
        "theta_b3": 1.0,
        "phi1": 1.5 + 1.5 * q,
        "phi2": 1.0 - q,
        "phi3": 0.0,
    }
    return np.array([vals[name] for name in PARAM_NAMES], dtype=float)


def t2_loose_vector(q, eta):
    vec = vector_for("T2")
    vec[PARAM_NAMES.index("phi2")] = 1.0 - q
    vec[PARAM_NAMES.index("theta_b2")] = q + eta
    return vec


def continuation_centers(channel, xvals, varied_name, center_name, search_half_width, samples):
    base = vector_for(channel)
    vi = PARAM_NAMES.index(varied_name)
    ci = PARAM_NAMES.index(center_name)
    centers = []
    offsets = np.linspace(-search_half_width, search_half_width, samples)
    for x in xvals:
        best_val = base[ci]
        best_p = -1.0
        for off in offsets:
            vec = base.copy()
            vec[vi] = base[vi] + x
            vec[ci] = base[ci] + off
            p = target_probability(channel, vec)
            if p > best_p:
                best_p = p
                best_val = vec[ci]
        centers.append(best_val)
    return np.array(centers)


def detuning_centers_for_t3(xvals, search_half_width=0.004, samples=121):
    base = vector_for("T3")
    phi2_i = PARAM_NAMES.index("phi2")
    delta_i = PARAM_NAMES.index("Delta")
    offsets = np.linspace(-search_half_width, search_half_width, samples)
    centers = []
    for x in xvals:
        best_delta = base[delta_i]
        best_p = -1.0
        for off in offsets:
            vec = base.copy()
            vec[phi2_i] = base[phi2_i] + x
            vec[delta_i] = base[delta_i] + off
            p = target_probability("T3", vec)
            if p > best_p:
                best_p = p
                best_delta = vec[delta_i]
        centers.append(best_delta)
    return np.array(centers)


def grid_from_rule(channel, xvals, yvals, rule):
    grid = np.empty((len(yvals), len(xvals)), dtype=float)
    for iy, y in enumerate(yvals):
        for ix, x in enumerate(xvals):
            grid[iy, ix] = target_probability(channel, rule(ix, x, y))
    return grid


def metrics(grid):
    return {
        "max": float(np.nanmax(grid)),
        "A99": float(np.nanmean(grid > 0.99)),
        "A999": float(np.nanmean(grid > 0.999)),
    }


def fmt_area(x):
    if x == 0:
        return "0"
    if x < 0.005:
        return f"{x:.1e}"
    return f"{x:.3f}"


def add_contours(ax, xs, ys, grid):
    levels = [level for level in (0.9, 0.99, 0.999) if np.nanmin(grid) <= level <= np.nanmax(grid)]
    if levels:
        colors = {0.9: "#ffffff", 0.99: "#ffd166", 0.999: "#7dd3fc"}
        widths = {0.9: 0.45, 0.99: 0.55, 0.999: 0.7}
        ax.contour(xs, ys, grid, levels=levels,
                   colors=[colors[l] for l in levels],
                   linewidths=[widths[l] for l in levels])


def annotate(ax, m):
    txt = rf"$P_\max={m['max']:.3f}$" + "\n" + rf"$A_{{.99}}={fmt_area(m['A99'])}$" + "\n" + rf"$A_{{.999}}={fmt_area(m['A999'])}$"
    ax.text(0.97, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.5, color="white",
            bbox=dict(boxstyle="round,pad=0.22", facecolor="black", alpha=0.55, edgecolor="none"))


def make_panels(n=170):
    panels = []

    # R1: direct-oracle local continuation in theta_a1 as phi2 is displaced.
    r1_base = vector_for("R1")
    x_r1 = np.linspace(-0.08, 0.08, n)
    y_r1 = np.linspace(-0.08, 0.08, n)
    phi2_i = PARAM_NAMES.index("phi2")
    theta_a1_i = PARAM_NAMES.index("theta_a1")
    r1_centers = continuation_centers("R1", x_r1, "phi2", "theta_a1", 0.12, 121)
    r1_naive = grid_from_rule(
        "R1", x_r1, y_r1,
        lambda ix, x, y: np.array([
            v if j not in (phi2_i, theta_a1_i) else (r1_base[phi2_i] + x if j == phi2_i else r1_base[theta_a1_i] + y)
            for j, v in enumerate(r1_base)
        ], dtype=float),
    )
    r1_guided = grid_from_rule(
        "R1", x_r1, y_r1,
        lambda ix, x, y: np.array([
            v if j not in (phi2_i, theta_a1_i) else (r1_base[phi2_i] + x if j == phi2_i else r1_centers[ix] + y)
            for j, v in enumerate(r1_base)
        ], dtype=float),
    )
    panels.append(("R1", x_r1, y_r1, r1_naive, r1_guided,
                   r"$R_1$: composite route",
                   r"$\delta\varphi_2/\pi$", r"$\delta\theta_{a1}/\pi$",
                   "numerical direct-oracle continuation of theta_a1"))

    # T2: analytic C2(q) branch versus a local q, eta scan.
    t2_base = vector_for("T2")
    q0 = 1.0 - t2_base[PARAM_NAMES.index("phi2")]
    eta0 = t2_base[PARAM_NAMES.index("theta_b2")] - q0
    x_t2 = np.linspace(0.0005, 0.06, n)
    y_t2 = np.linspace(-0.02, 0.02, n)
    x_t2[np.argmin(np.abs(x_t2 - q0))] = q0
    y_t2[np.argmin(np.abs(y_t2 - eta0))] = eta0
    t2_naive = grid_from_rule("T2", x_t2, y_t2, lambda ix, q, eta: t2_loose_vector(q, eta))
    t2_guided = grid_from_rule("T2", x_t2, y_t2, lambda ix, q, eta: c2_vector(q, eta))
    panels.append(("T2", x_t2, y_t2, t2_naive, t2_guided,
                   r"$T_2$: shared-waveguide route",
                   r"$q$", r"$\eta/\pi$",
                   "analytic C2(q) branch"))

    # T3: detuning-guided continuation around the near-resonant optimized point.
    t3_base = vector_for("T3")
    x_t3 = np.linspace(-0.006, 0.006, n)
    y_t3 = np.linspace(-0.006, 0.006, n)
    phi2_i = PARAM_NAMES.index("phi2")
    theta_b2_i = PARAM_NAMES.index("theta_b2")
    delta_i = PARAM_NAMES.index("Delta")
    delta_centers = detuning_centers_for_t3(x_t3)
    t3_naive = grid_from_rule(
        "T3", x_t3, y_t3,
        lambda ix, x, y: np.array([
            v if j not in (phi2_i, theta_b2_i) else (t3_base[phi2_i] + x if j == phi2_i else t3_base[theta_b2_i] + y)
            for j, v in enumerate(t3_base)
        ], dtype=float),
    )
    t3_guided = grid_from_rule(
        "T3", x_t3, y_t3,
        lambda ix, x, y: np.array([
            (delta_centers[ix] if j == delta_i else
             t3_base[phi2_i] + x if j == phi2_i else
             t3_base[theta_b2_i] + y if j == theta_b2_i else v)
            for j, v in enumerate(t3_base)
        ], dtype=float),
    )
    panels.append(("T3", x_t3, y_t3, t3_naive, t3_guided,
                   r"$T_3$: near-resonant route",
                   r"$\delta\varphi_2/\pi$", r"$\delta\theta_{b2}/\pi$",
                   "numerical detuning-guided direct-oracle continuation"))
    return panels


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(REPORT_DIR, exist_ok=True)
    panels = make_panels()

    fig, axes = plt.subplots(3, 2, figsize=(7.1, 7.2), constrained_layout=True)
    last_im = None
    report_rows = []
    for row, (channel, xs, ys, naive, guided, row_title, xlabel, ylabel, method) in enumerate(panels):
        for col, (grid, col_title) in enumerate([(naive, "unconstrained local scan"), (guided, "mechanism-guided scan")]):
            ax = axes[row, col]
            last_im = ax.pcolormesh(xs, ys, grid, shading="auto", cmap=NIGHTFIRE, vmin=0, vmax=1, rasterized=True)
            add_contours(ax, xs, ys, grid)
            m = metrics(grid)
            annotate(ax, m)
            ax.set_title(f"{row_title}: {col_title}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            report_rows.append((channel, col_title, method if col else "local two-parameter scan around optimized point", m))

    cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(), shrink=0.94, pad=0.012)
    cbar.set_label("target probability")
    for basename in ("figure3_mechanism_comparison", "fig3_correlated_branch"):
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(OUT_DIR, f"{basename}.{ext}"))
    plt.close(fig)

    report_path = os.path.join(REPORT_DIR, "figure3_mechanism_comparison_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Figure 3 Mechanism Comparison Report\n\n")
        f.write("This figure compares representative mechanisms for `R1`, `T2`, and `T3`. It intentionally omits `R2` and `R3` from the main comparison because `T2/R2` and `T3/R3` are direction-related pairs under the fixed-backbone symmetry.\n\n")
        f.write("| Row | Panel | Basis | Max target probability | A(P>0.99) | A(P>0.999) |\n")
        f.write("| --- | --- | --- | --- | --- | --- |\n")
        for channel, panel, method, m in report_rows:
            f.write(f"| {channel} | {panel} | {method} | {m['max']:.6f} | {m['A99']:.6f} | {m['A999']:.6f} |\n")
        f.write("\n## Scan definitions\n\n")
        f.write("- `R1` row: scans `delta phi2/pi` and `delta theta_a1/pi` around the optimized `R1` point. The guided panel uses numerical direct-oracle continuation: for each `delta phi2/pi`, `theta_a1` is recentered by a one-dimensional direct-oracle maximization, then the transverse offset is scanned.\n")
        f.write("  Held fixed except for the scanned/recentered parameters: `Gamma`, `Delta`, `theta_a2`, `theta_b2`, `theta_b3`, `phi1`, and `phi3` at the optimized `R1` point. No analytic branch is imposed.\n")
        f.write("- `T2` row: scans `(q, eta/pi)`. The unconstrained panel varies only local `(q, eta)` around the optimized point; the guided panel imposes the analytic `C2(q)` branch. This is the only row using a closed analytic branch.\n")
        f.write("  Guided constraints: `phi2/pi=1-q`, `phi3/pi=0`, `theta_b3/pi=1`, `phi1/pi=3/2+3q/2`, `theta_a1/pi=1/2+3q/2`, `Delta/Gamma=2 sin(pi q)`, and `theta_b2/pi=q+eta`; `Gamma=1` and `theta_a2/pi=1` are held fixed.\n")
        f.write("- `T3` row: scans `delta phi2/pi` and `delta theta_b2/pi` around the optimized near-resonant `T3` point. The guided panel uses numerical direct-oracle continuation of `Delta`: for each `delta phi2/pi`, `Delta` is recentered by a one-dimensional maximization, then the transverse `theta_b2` offset is scanned.\n")
        f.write("  Held fixed except for the scanned/recentered parameters: `Gamma`, `theta_a1`, `theta_a2`, `theta_b3`, `phi1`, and `phi3` at the optimized `T3` point. No analytic branch is imposed.\n")
        f.write("\nAll probabilities are computed with the direct 14x14 oracle for WG1-L incidence; no deprecated one-factor formulas for `t2`, `t3`, or `r3` are used.\n")


if __name__ == "__main__":
    main()
