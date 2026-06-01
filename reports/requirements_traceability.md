# Requirements Traceability

| Requirement | Status | Evidence |
|---|---|---|
| Abstract terminology simplified and undefined terms removed | Done | `paper/main.tex`, Abstract: uses direct descriptions of routing, eta-form, and time-reversal partner relation without early unexplained jargon. |
| Appendix and supplementary material clearly distinguished | Done | `paper/main.tex`: no `Supplemental Material` claims; same-PDF derivations are called Appendix and code is called online repository. |
| Port naming simplified and not unnecessarily indirect | Done | `paper/main.tex`, Eq. (1) and surrounding text define `(WG1-L, WG1-R, WG2-L, WG2-R, WG3-L, WG3-R)` directly. |
| Eta notation explained clearly, including physical meaning and why it is introduced | Done | `paper/main.tex`, Hamiltonian section defines `eta_m=g f_m/v_g`; edge-incidence and middle-incidence sections state that eta amplitudes are atomic excitation amplitudes used to express outgoing channels compactly. |
| Phase-factor notation is consistent | Done | `paper/main.tex`, phase-factor definitions around Eqs. for `P_{j,m}` and trigonometric forms; scripts use the same `PARAM_NAMES` in `scripts/config.py`. |
| `N_{t1}` described as numerator abbreviation, not a physical factorization | Done | `paper/main.tex`, WG1-L discussion identifies `N_{t_1}` as the full numerator for `t_1`; Table I treats it as a zero condition, not as an elementary phase factor. |
| Fig. 2 follows representative-channel requirement: `T1`, `R1`, `T2`, `T3` only | Done | `paper/main.tex`, Fig. 2 caption and text; generated `figures/fig2_representative_maps.pdf`; script `scripts/regenerate_fig2.py` uses `PLOTS_MAIN = [T1,R1,T2,T3]`. |
| Delta-sensitivity discussion and figure are present | Done | `paper/main.tex`, `Detuning sensitivity of near-resonant routes`; `figures/fig4_detuning_sensitivity.pdf`; script `scripts/regenerate_fig4.py`. |
| Small-atom comparison is present and connected to giant-atom advantage | Done | `paper/main.tex`, Fig. 1 bottom panel and `Small-atom reference` subsection explain absence of independent directional phase zeros. |
| Middle-incidence equations explicit enough for back-substitution checks | Done | `paper/main.tex`, middle-incidence amplitude formulas and Appendix 2x2 reduction; `symbolic/middle_incidence_2x2_system.tex`; `audit/audit_WG2L_eta_form.py`. |
| Repository contains code and formulas needed to verify the 6x6 matrix | Done | `src/direct_solver_6x6.py`, `src/validate_formulas.py`, `src/symmetry_checks.py`, `audit/run_formula_audit.py`, `symbolic/eta_expressions.tex`. |
| Section IV.B and IV.D explain method traditionally: ansatz, jump equations, atomic equations, 14x14 system, reduction to 2x2, back-substitution verification | Done | `paper/main.tex`: edge and middle scattering formulations, Appendix explicit equations, 2x2 middle-incidence system, Verification subsection; scripts `audit/audit_WG2L_eta_form.py` and `audit/run_formula_audit.py`. |
| `C_2(q)` introduced only where needed and not used prematurely | Done | `paper/main.tex`, `C_2(q)` appears in the dedicated correlated-branch subsection after representative Fig. 2 is introduced; Fig. 2 caption refers only to a branch neighborhood introduced below. |
| Main Fig. 2 does not contain `R2` | Done | `scripts/regenerate_fig2.py`, `PLOTS_MAIN`; rendered `figures/fig2_representative_maps.pdf`. |
| Main Fig. 2 does not contain `R3` | Done | `scripts/regenerate_fig2.py`, `PLOTS_MAIN`; rendered `figures/fig2_representative_maps.pdf`. |
| Full six-target data retained outside main figure | Done | `figures/fig2_all6_supplementary.pdf`, `data/optimized_params.json`, Table III in `paper/main.tex`. |

No item in this report is marked Done without an implementation file, manuscript location, figure, table, or script.
