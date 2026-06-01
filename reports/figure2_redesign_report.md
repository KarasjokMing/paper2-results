# Figure 2 Redesign Report

## Main figure

The revised main Fig. 2 is `manuscript/figures/fig_phase_routing_maps.pdf`, also exported as `figure2_representative.pdf`.
It contains exactly four representative panels:

| Panel | Role | Coordinates |
| --- | --- | --- |
| `T1` | Direct self-transmission / sufficient entrance phase-zero route | fixed `phi1/pi=0.25`, fixed `theta_a1/pi=0.75`, scanned `phi2/pi` and `theta_b2/pi` over `[0,2] x [0,2]` |
| `R1` | Same-waveguide reflection / composite route | scanned `phi2/pi` and `theta_a1/pi` |
| `T2` | Representative shared-waveguide route | scanned `phi1/pi` and `theta_a1/pi` in an extended branch window |
| `T3` | Representative far-waveguide near-resonant route | scanned `phi2/pi` and `theta_b2/pi` in a magnified near-resonant window |

The mandatory `T1` redesign is implemented in `03_repro_scripts/fine_search_core/make_global_phase_maps.py` and `_github_upload/scripts/regenerate_fig2.py` through the `t1_sufficient_map` function.
It imposes
\[
P_{1+}=1+\exp[i(\phi_1+\theta_{a1})]=0
\]
by fixing `phi1/pi=0.25` and `theta_a1/pi=0.75`, then scans the unrelated downstream phases `phi2/pi` and `theta_b2/pi`.
This demonstrates sufficiency: once the entrance interference condition is imposed, \(T_1\) remains near unity over a broad transverse phase plane except possible singular denominator points.
It does not claim that \(P_{1+}=0\) is necessary for high \(T_1\).

## Supplementary/repository all-target figure

The complete six-panel record is retained outside the main text as `figure2_all6_supplementary.pdf`.
It contains `T1`, `R1`, `T2`, `R2`, `T3`, and `R3`, preserving the full optimized target display for repository inspection.
The main manuscript omits `R2` and `R3` because `T2/R2` and `T3/R3` are direction-related pairs under the fixed-backbone symmetry.

## Generated files

- `figure2_representative.pdf`
- `figure2_all6_supplementary.pdf`
- `t1_sufficient_condition_scan.pdf`
- `manuscript/figures/fig_phase_routing_maps.pdf`
- `manuscript/figures/fig_phase_routing_maps_all6.pdf`
- `manuscript/figures/t1_sufficient_condition_scan.pdf`

## Reproduction

Project manuscript figure:

```bash
python 03_repro_scripts/fine_search_core/make_global_phase_maps.py
```

Repository figure:

```bash
cd _github_upload
python scripts/regenerate_fig2.py
```
