# Middle-Incidence Scope Report

The revised manuscript uses the minimal-scope option for middle-waveguide incidence.

## Scope Chosen

The paper keeps the exact transparent middle-incidence point \(\tilde t_2=1\) and does not add a new six-target optimization table for
\[
\tilde r_1,\tilde t_1,\tilde r_2,\tilde t_2,\tilde r_3,\tilde t_3.
\]
This choice is stated explicitly in `manuscript/main.tex`, where the text says that a systematic search for the remaining five middle-incidence target outputs is not pursued in the present manuscript and is left for repository-level numerical exploration.

## Equation Support

The middle-incidence derivation is still fully checkable:

- `manuscript/main.tex` contains the WG2 ansatz, jump conditions, \(\eta\)-form amplitudes, \(2\times2\) atomic reduction, and back-substitution verification.
- `middle_incidence_equations_full.tex` records the full 14-equation middle-incidence system in standalone form.
- `_github_upload/audit/audit_WG2L_eta_form.py` verifies the WG2 eta-form against the direct \(14\times14\) solver.
- `_github_upload/reports/validation_report.md` and `_github_upload/reports/formula_audit_report.md` summarize validation results.

## Caveat

The manuscript does not claim complete analytical routing classifications for the five nontransparent middle-incidence targets.
Those routes are treated as future numerical exploration using the same direct-oracle framework.
