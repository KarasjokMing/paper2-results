"""Reproduce the detuning-sensitivity figure used as Fig. 4 in the repository package."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import CHANNELS, PARAM_NAMES
from physics import compute_numpy


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "paper" / "figures"
PARAM_PATH = Path(__file__).with_name("optimal_params.json")


def probabilities(channel: str, delta_values: np.ndarray) -> np.ndarray:
    with PARAM_PATH.open("r", encoding="utf-8") as f:
        opt = json.load(f)
    vec = np.array([opt[channel][name] for name in PARAM_NAMES], dtype=float)
    idx = CHANNELS.index(channel)
    out = np.empty_like(delta_values, dtype=float)
    for i, delta in enumerate(delta_values):
        vec[PARAM_NAMES.index("Delta")] = delta
        p = compute_numpy(vec[0], vec[1], *vec[2:])
        out[i] = p[idx] if np.all(np.isfinite(p)) else np.nan
    return out


def best_at_fixed_delta(channel: str, fixed_delta: float | None) -> float:
    with PARAM_PATH.open("r", encoding="utf-8") as f:
        opt = json.load(f)
    vec = np.array([opt[channel][name] for name in PARAM_NAMES], dtype=float)
    if fixed_delta is not None:
        vec[PARAM_NAMES.index("Delta")] = fixed_delta
    p = compute_numpy(vec[0], vec[1], *vec[2:])
    return float(p[CHANNELS.index(channel)])


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with PARAM_PATH.open("r", encoding="utf-8") as f:
        opt = json.load(f)
    delta_star = float(opt["T3"]["Delta"])
    delta = np.linspace(delta_star - 0.012, delta_star + 0.012, 1200)
    t3 = probabilities("T3", delta)

    channels = ["T1", "R1", "T2", "R2", "T3", "R3"]
    free = [best_at_fixed_delta(ch, None) for ch in channels]
    resonant = [best_at_fixed_delta(ch, 0.0) for ch in channels]

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax1.plot(delta, t3, color="#2166ac", lw=1.0)
    ax1.axvline(delta_star, color="#b2182b", ls="--", lw=0.8)
    ax1.axhline(0.999, color="#d6604d", ls=":", lw=0.8)
    ax1.set_xlabel(r"$\Delta/\Gamma$")
    ax1.set_ylabel(r"$T_3$")
    ax1.set_title(r"$T_3$ detuning window")
    ax1.set_ylim(0, 1.04)

    x = np.arange(len(channels))
    width = 0.36
    ax2.bar(x - width / 2, free, width, label=r"optimized $\Delta$", color="#2166ac")
    ax2.bar(x + width / 2, resonant, width, label=r"$\Delta=0$", color="#d6604d")
    ax2.set_xticks(x, channels)
    ax2.set_ylim(0, 1.04)
    ax2.set_ylabel("target probability")
    ax2.set_title("strict-resonance comparison")
    ax2.legend(frameon=False, fontsize=7)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig6_detuning_degeneracy.pdf")
    fig.savefig(OUT_DIR / "fig6_detuning_degeneracy.png")
    print(f"Wrote {OUT_DIR / 'fig6_detuning_degeneracy.pdf'}")


if __name__ == "__main__":
    main()
