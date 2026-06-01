# WG2-L Eta-Form Audit

Source script: `audit/audit_WG2L_eta_form.py`.

Summary:
- The WG2-L eta-form was verified by direct back-substitution into the 14 Schrodinger equations.
- The analytic amplitudes were compared against the numerical 14x14 solver on random parameter points under the fixed backbone `Gamma=1`, `theta_a2=pi`.
- Representative maximum errors are at double-precision roundoff scale. In the included audit run, the 14-equation back-substitution residual was about `2.1e-15`, the analytic-vs-solver amplitude difference was about `5.4e-15`, and the probability-conservation deviation was about `1.3e-14`.

Command:

```bash
python audit/audit_WG2L_eta_form.py
```
