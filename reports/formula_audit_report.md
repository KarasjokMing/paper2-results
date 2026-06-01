# Formula Audit Report

Commands run from repository root:

```bash
python src/validate_formulas.py
python src/symmetry_checks.py
python audit/run_formula_audit.py
```

Summary:

| Check | Result | Maximum residual/error |
|---|---:|---:|
| Probability conservation for all six direct-solver columns | PASS | `2.66e-15` |
| 6x6 unitarity `||S^\dagger S-I||_\infty` | PASS | `2.66e-15` |
| Time-reversal partner relation | PASS | `3.55e-15` |
| Zero-phase reciprocity | PASS | `1.89e-15` |
| WG1-L eta-form vs direct solver | PASS | `1.67e-15` |
| WG2-L eta-form back-substitution into 14 equations | PASS | `2.07e-15` |
| WG2-L eta-form vs 14x14 solver | PASS | `5.39e-15` |
| All six incident ports eta-form vs direct solver | PASS | `1.96e-15` |

Implementation files:

- Direct matrix construction: `src/direct_solver_6x6.py`.
- WG1-L eta-form audit: `audit/audit_WG1L_eta_form.py`.
- WG2-L eta-form audit and explicit back-substitution: `audit/audit_WG2L_eta_form.py`.
- All-port eta-form audit: `audit/audit_all_ports_eta_form.py` and `src/symmetry_checks.py`.

All values are consistent with double-precision roundoff.
