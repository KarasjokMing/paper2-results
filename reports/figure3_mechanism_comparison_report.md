# Figure 3 Mechanism Comparison Report

This figure compares representative mechanisms for `R1`, `T2`, and `T3`. It intentionally omits `R2` and `R3` from the main comparison because `T2/R2` and `T3/R3` are direction-related pairs under the fixed-backbone symmetry.

| Row | Panel | Basis | Max target probability | A(P>0.99) | A(P>0.999) |
| --- | --- | --- | --- | --- | --- |
| R1 | unconstrained local scan | local two-parameter scan around optimized point | 0.999998 | 0.423426 | 0.129412 |
| R1 | mechanism-guided scan | numerical direct-oracle continuation of theta_a1 | 0.999998 | 0.423460 | 0.129412 |
| T2 | unconstrained local scan | local two-parameter scan around optimized point | 1.000000 | 0.000069 | 0.000035 |
| T2 | mechanism-guided scan | analytic C2(q) branch | 1.000000 | 0.213356 | 0.059135 |
| T3 | unconstrained local scan | local two-parameter scan around optimized point | 1.000000 | 0.086332 | 0.046644 |
| T3 | mechanism-guided scan | numerical detuning-guided direct-oracle continuation | 1.000000 | 0.086332 | 0.046644 |

## Scan definitions

- `R1` row: scans `delta phi2/pi` and `delta theta_a1/pi` around the optimized `R1` point. The guided panel uses numerical direct-oracle continuation: for each `delta phi2/pi`, `theta_a1` is recentered by a one-dimensional direct-oracle maximization, then the transverse offset is scanned.
  Held fixed except for the scanned/recentered parameters: `Gamma`, `Delta`, `theta_a2`, `theta_b2`, `theta_b3`, `phi1`, and `phi3` at the optimized `R1` point. No analytic branch is imposed.
- `T2` row: scans `(q, eta/pi)`. The unconstrained panel varies only local `(q, eta)` around the optimized point; the guided panel imposes the analytic `C2(q)` branch. This is the only row using a closed analytic branch.
  Guided constraints: `phi2/pi=1-q`, `phi3/pi=0`, `theta_b3/pi=1`, `phi1/pi=3/2+3q/2`, `theta_a1/pi=1/2+3q/2`, `Delta/Gamma=2 sin(pi q)`, and `theta_b2/pi=q+eta`; `Gamma=1` and `theta_a2/pi=1` are held fixed.
- `T3` row: scans `delta phi2/pi` and `delta theta_b2/pi` around the optimized near-resonant `T3` point. The guided panel uses numerical direct-oracle continuation of `Delta`: for each `delta phi2/pi`, `Delta` is recentered by a one-dimensional maximization, then the transverse `theta_b2` offset is scanned.
  Held fixed except for the scanned/recentered parameters: `Gamma`, `theta_a1`, `theta_a2`, `theta_b3`, `phi1`, and `phi3` at the optimized `T3` point. No analytic branch is imposed.

All probabilities are computed with the direct 14x14 oracle for WG1-L incidence; no deprecated one-factor formulas for `t2`, `t3`, or `r3` are used.
