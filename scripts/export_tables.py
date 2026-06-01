"""Export compact numerical tables used by the manuscript."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARAMS = ["Gamma", "Delta", "theta_a1", "theta_a2", "theta_b2", "theta_b3", "phi1", "phi2", "phi3"]


def main() -> None:
    with (ROOT / "data" / "optimized_params.json").open("r", encoding="utf-8") as f:
        data = json.load(f)

    out = ROOT / "reports" / "optimized_parameter_table.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("# Optimized Parameter Table\n\n")
        f.write("| Target | " + " | ".join(PARAMS) + " |\n")
        f.write("|---|" + "|".join(["---"] * len(PARAMS)) + "|\n")
        for target, vals in data.items():
            row = [target] + [f"{float(vals[p]):.10g}" for p in PARAMS]
            f.write("| " + " | ".join(row) + " |\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
