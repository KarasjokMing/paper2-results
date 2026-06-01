# Giant-Atom Six-Port Scattering

This repository contains the reproducibility package for the manuscript on single-photon routing in a two-giant-atom, three-waveguide node. It includes the checked manuscript, direct 14x14/6x6 solvers, eta-form audits, figure-regeneration scripts, symbolic derivation files, optimized parameters, and traceability reports.

## Port Convention

The six-port scattering matrix is ordered as

```text
(WG1-L, WG1-R, WG2-L, WG2-R, WG3-L, WG3-R)
```

For WG1-L incidence, the manuscript reports the output vector as

```text
(R1, T1, R2, T2, R3, T3)
```

The direct solver constructs

```text
S_dir = P A^{-1} B
```

where `A` is the 14x14 coefficient matrix, `B` stores the six incident boundary conditions, and `P` projects the internal solution onto the six output ports.

## Main Files

- `paper/main.tex`, `paper/main.pdf`: checked manuscript source and PDF.
- `src/direct_oracle.py`: direct 14x14 solver and 6x6 scattering matrix construction.
- `src/direct_solver_6x6.py`: same implementation retained under its development filename.
- `src/formulas_eta.py`: numerical eta-form probability evaluator used by the figure scripts.
- `src/validate_formulas.py`: probability conservation, unitarity, and time-reversal partner checks.
- `src/symmetry_checks.py`: all-port eta-form comparison against the direct solver.
- `audit/run_formula_audit.py`: complete formula-audit runner, including WG1-L, WG2-L, all-port, and direct-oracle checks.

## Figure Scripts

Run from the repository root:

```bash
python scripts/regenerate_fig2.py
python scripts/regenerate_fig3.py
python scripts/regenerate_fig4.py
```

Outputs:

- `figures/figure2_representative.pdf`: main-text Fig. 2 with exactly `T1`, `R1`, `T2`, `T3`; the `T1` panel fixes `P_{1+}=0` and scans unrelated phases.
- `figures/figure2_all6_supplementary.pdf`: all-six repository map containing `T1`, `R1`, `T2`, `R2`, `T3`, `R3`.
- `figures/t1_sufficient_condition_scan.pdf`: standalone `T1` robustness scan with `P_{1+}=0` imposed.
- `figures/fig3_correlated_branch.pdf`: correlated `C_2(q)` branch verification.
- `figures/fig4_detuning_sensitivity.pdf`: near-resonant detuning-sensitivity comparison.

## Validation Commands

Run from the repository root:

```bash
python src/validate_formulas.py
python src/symmetry_checks.py
python audit/run_formula_audit.py
```

Expected residuals:

- probability conservation for all six columns: below about `3e-15`;
- unitarity of the 6x6 scattering matrix: below about `3e-15`;
- time-reversal partner relation: below about `4e-15`;
- eta-form vs direct 14x14 solver: below about `3e-15` for all ports;
- explicit 14-equation back-substitution residuals: about `2e-15` for WG2-L.

## Manuscript Mapping

- Eq. `S_dir = P A^{-1} B`: `src/direct_oracle.py`.
- WG1-L eta-form checks: `audit/audit_WG1L_eta_form.py`.
- WG2-L eta-form and back-substitution checks: `audit/audit_WG2L_eta_form.py`.
- All six incident ports: `src/symmetry_checks.py` and `audit/audit_all_ports_eta_form.py`.
- Fig. 2: `scripts/regenerate_fig2.py`.
- Fig. 3: `scripts/regenerate_fig3.py`.
- Fig. 4 detuning sensitivity: `scripts/regenerate_fig4.py`.
- Optimized parameter table: `scripts/export_tables.py`.

## Reports

- `reports/validation_report.md`: validation summary and acceptance checks.
- `reports/formula_audit_report.md`: formula-audit summary.
- `reports/requirements_traceability.md`: requirement-by-requirement traceability.
- `reports/author_comments_traceability.md`: traceability for the 13 author comments addressed in the manuscript revision.
- `reports/figure2_redesign_report.md`: parameters and design rationale for revised Fig. 2.
- `reports/middle_incidence_scope_report.md`: scope statement for the middle-incidence treatment.

## Deprecated Code

`python/deprecated/oracle_compact.py` is retained only for provenance. It constructs the old symmetry-completed matrix using `S=S^T` and waveguide reversal. It is not the exact 6x6 scattering matrix and must not be used as the truth source. Use `src/direct_oracle.py` or `src/direct_solver_6x6.py` instead.
