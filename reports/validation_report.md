# Validation Report

Generated after running:

```bash
python src/validate_formulas.py
python src/symmetry_checks.py
python audit/run_formula_audit.py
```

Required validation items:

| Required item | Status | Evidence |
|---|---|---|
| Probability conservation for all six columns | PASS | `src/validate_formulas.py`: max column-conservation deviation `2.66e-15`. |
| Unitarity of the 6x6 matrix | PASS | `src/validate_formulas.py`: max `||S^\dagger S-I||_\infty = 2.66e-15`. |
| Eta-form formulas vs direct 14x14 solver | PASS | `src/symmetry_checks.py`: all six ports pass, max error `1.96e-15`; `audit/audit_WG1L_eta_form.py`: max WG1-L error `1.67e-15`; `audit/audit_WG2L_eta_form.py`: WG2-L analytic-vs-solver error `5.39e-15`. |
| Time-reversal partner relation | PASS | `src/validate_formulas.py`: max TR-partner deviation `3.55e-15`. |
| Back-substitution residuals for explicit equations | PASS | `audit/audit_WG2L_eta_form.py`: 14-equation back-substitution residual `2.07e-15`; manuscript Appendix gives the explicit 14-equation system. |

Figure reproduction:

| Figure | Command | Status |
|---|---|---|
| Main Fig. 2 representative maps | `python scripts/regenerate_fig2.py` | PASS; writes `figures/fig2_representative_maps.pdf` with `T1`, `R1`, `T2`, `T3` only. |
| Repository all-six map | `python scripts/regenerate_fig2.py` | PASS; writes `figures/fig2_all6_supplementary.pdf`. |
| Fig. 3 correlated branch | `python scripts/regenerate_fig3.py` | PASS; writes `figures/fig3_correlated_branch.pdf`. |
| Fig. 4 detuning sensitivity | `python scripts/regenerate_fig4.py` | PASS; writes `figures/fig4_detuning_sensitivity.pdf`. |

Conclusion: the repository contains the code, optimized data, formulas, and scripts needed to reproduce the manuscript validation and the requested representative Fig. 2.
