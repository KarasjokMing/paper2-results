# Author Comments Traceability

| Author comment | Status: Done/Partial/Not done | Exact manuscript change | Figure/script/repository evidence | Remaining caveat |
| --- | --- | --- | --- | --- |
| 1. Simplify abstract | Done | `manuscript/main.tex:40-46` rewrites the abstract around the model, geometry, exact solution, full \(6\times6\) matrix, phase-controlled routing, small-atom contrast, and plain time-reversal relation. | Compiled in `main_author_comments_revised.pdf`. | None. |
| 2. Distinguish Appendix, Supplemental Material, and GitHub | Done | `manuscript/main.tex` uses Appendix for in-PDF derivations and online repository for code/data; no separate Supplemental Material is cited. The verification-code statement names `https://github.com/KarasjokMing/paper2-results`, branch `codex/repro-package`, commit `63a9204`. | Search performed for `Supplemental Material`, `supplementary`, `repository`, and `machine-readable`; remaining repository statements are online-repository statements. | None. |
| 3. Direct port convention | Done | `manuscript/main.tex:101-110` defines `(1,2,3,4,5,6)` and immediately states \(S_{ij}\) is output port \(i\) for input port \(j\); `s^{(1)}` is given as `(r1,t1,r2,t2,r3,t3)^T`. | Compiled in `main_author_comments_revised.pdf`. | None. |
| 4. Explain eta-form and reconcile \(r_2\) with eta-form | Done | `manuscript/main.tex:208-213` explains \(\eta_a^{(j)},\eta_b^{(j)}\), incident-port superscripts, source dependence, and separation of atomic response from radiation factors; `manuscript/main.tex:255-257` states \(r_2=-\Gamma P_{1+}A_r/M\) is verified and the simple \(t_2,t_3,r_3\) eliminated forms failed audit. | `_github_upload/audit/audit_WG1L_eta_form.py`; `_github_upload/reports/formula_audit_report.md`; `_github_upload/outputs/reports/WG1L_eta_form_audit.md`. | None. |
| 5. Consistent phase-factor notation | Done | `manuscript/main.tex:225-234` defines \(P_{j,m}^{\rm src}\), \(P_{j,m}^{L}\), and \(P_{j,m}^{R}\), explains \(P^R=(P^{\rm src})^*\), and defines \(P_{1+},P_{1-},P_{3+},P_{3-}\) as shorthand. | Search removed unexplained `P_{1+}^*`, `P_{2+}^{a*}`, and `P_{3+}^*` from the main amplitude discussion. | Some shorthand remains where already defined and physically useful. |
| 6. Describe \(N_{t1}\) correctly and keep factorized section concise | Done | `manuscript/main.tex:250-257` states \(N_{t_1}\) is only the full numerator of \(t_1\), with explicit expression in the Appendix, and lists which amplitudes factorize versus remain in verified eta-form. | Appendix expression retained at `manuscript/main.tex:702-718`. | None. |
| 7. Old branch/phase scatter plot / Fig. 3 mechanism comparison | Done | No disconnected old branch/scatter plot was added. The former T2/R2-only Fig. 3 text was replaced by a mechanism-guided comparison of \(R_1,T_2,T_3\) in `manuscript/main.tex`, subsection "Mechanism-guided basin comparison". | `manuscript/figures/fig_constructive_branch_verification.pdf` now contains three rows \(R_1,T_2,T_3\); `_github_upload/scripts/regenerate_fig3_mechanism_comparison.py`; `_github_upload/figures/figure3_mechanism_comparison.pdf`; `_github_upload/reports/figure3_mechanism_comparison_report.md`. | The \(R_1\) and \(T_3\) guided scans are explicitly numerical direct-oracle continuations, not analytic branches. |
| 8. Restore detuning-sensitivity discussion and figure | Done | `manuscript/main.tex:429-445` discusses strict \(\Delta=0\) versus small nonzero \(\Delta\), keeps the detuning figure, and avoids a universal \(25/36\) claim. | `manuscript/figures/fig6_detuning_degeneracy.pdf`; `_github_upload/scripts/regenerate_fig4.py`; `_github_upload/figures/fig4_detuning_sensitivity.pdf`. | No \(25/36\) dashed line is present in the current manuscript figure, so no such caption is needed. |
| 9. Small-atom comparison and schematic | Done | `manuscript/main.tex:88-96` captions Fig. 1 as a two-panel giant/small-atom comparison; `manuscript/main.tex:448-480` explains absence of independent local phases and directional tuning in the small-atom reference. | `manuscript/figures/fig1_system_schematic_vector_cropped.pdf`; `manuscript/figures/fig_small_atom_schematic_vector_cropped.pdf`. | None. |
| 10. Middle incidence scope | Done | `manuscript/main.tex:667` states the remaining five middle-incidence target searches are not pursued and are left for repository-level numerical exploration. | `middle_incidence_scope_report.md`; `_github_upload/audit/audit_WG2L_eta_form.py`. | Full six-target middle-incidence optimization is intentionally not included. |
| 11. Fig. 2 T1 sufficient-condition scan | Done | `manuscript/main.tex:373-382` caption says the T1 panel imposes \(\Poneplus=0\), scans unrelated \((\varphi_2,\theta_{b2})\), and is a sufficiency test, not a necessity proof; `manuscript/main.tex:385-390` explains the plateau. | `03_repro_scripts/fine_search_core/make_global_phase_maps.py`; `_github_upload/scripts/regenerate_fig2.py`; `t1_sufficient_condition_scan.pdf`; `figure2_redesign_report.md`. | The panel is visually near-uniform by design because the sufficient condition makes \(T_1\) robust over the scanned plane. |
| 12. Fig. 2 representative channels only | Done | `manuscript/main.tex:358-363` explains why the main Fig. 2 shows only \(T_1,R_1,T_2,T_3\); caption at `manuscript/main.tex:372-382` refers to four mechanisms. | `figure2_representative.pdf` has exactly four panels; `figure2_all6_supplementary.pdf` retains all six for repository/supplementary inspection. | None. |
| 13. \(C_2(q)\) appears only in the mechanism-guided Fig. 3 discussion | Done | Section III.A discusses Fig. 2 qualitatively; \(\mathcal{C}_2(q)\) first appears in the Fig. 3 mechanism-guided subsection, where it is used only for the \(T_2\) row. | Search confirms no literal `\mathcal{C}_2(q)` before the Fig. 3 subsection. Fig. 3 includes \(R_1,T_2,T_3\), while \(R_2,R_3\) are omitted from the main mechanism comparison and retained in the table/repository. | None. |

## Final Acceptance Checklist

| Item | Status |
| --- | --- |
| Abstract simplified | Done |
| Appendix/Supplemental/GitHub terminology standardized | Done |
| Port convention direct | Done |
| Eta-form explanation added | Done |
| Phase notation consistent | Done |
| \(N_{t1}\) described as numerator abbreviation | Done |
| Small-atom comparison retained and schematic checked | Done |
| Delta-sensitivity discussion and figure retained | Done |
| Middle-incidence scope clarified | Done |
| Explicit middle-incidence equations/back-substitution support retained | Done |
| T1 sufficient-condition scan added | Done |
| Main Fig. 2 uses T1/R1/T2/T3 only | Done |
| \(C_2(q)\) introduced only in Fig. 3 branch discussion | Done |
| Fig. 3 includes R1, T2, and T3 rather than only T2/R2 | Done |
| GitHub repository populated and commit hash returned | Done: `63a9204` on branch `codex/repro-package` |
