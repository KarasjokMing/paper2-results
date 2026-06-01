# Giant-Atom Six-Port Scattering

This repository contains the reproducibility package for the manuscript on single-photon routing in a two-giant-atom, three-waveguide node. The package includes the checked manuscript, the exact direct-oracle solver, formula-audit scripts, figure-reproduction scripts, and symbolic expression summaries.

## Port Convention

The scattering matrix is ordered as

```text
(WG1-L, WG1-R, WG2-L, WG2-R, WG3-L, WG3-R)
```

For waveguide-1 left incidence, the output probabilities are reported as

```text
(R1, T1, R2, T2, R3, T3)
```

The exact direct solver constructs the full six-port matrix as

```text
S_dir = P A^{-1} B
```

where `A` is the 14x14 coefficient matrix, `B` contains the six incident boundary conditions, and `P` projects the internal solution onto the six output ports.

## Parameters

`Gamma` is the waveguide decay scale, `Delta` is the common atom-photon detuning, `phi_j` are propagation phases along waveguide `j`, and `theta_a1`, `theta_a2`, `theta_b2`, `theta_b3` are synthetic local coupling phase differences. The manuscript uses the fixed backbone `Gamma=1` and `theta_a2=pi` for the representative routing landscape.

## Validation Commands

Run from the repository root:

```bash
python audit/validate_direct_oracle.py
python audit/audit_WG1L_eta_form.py
python audit/audit_WG2L_eta_form.py
python audit/audit_all_ports_eta_form.py
python audit/run_formula_audit.py
```

Expected residuals are at double-precision roundoff scale: direct-oracle probability conservation and unitarity below approximately `3e-15`, eta-form comparisons below approximately `3e-15`, and 14-equation back-substitution residuals below approximately `2e-15`.

## Figure Commands

Run from the repository root:

```bash
python figures/reproduce_fig2.py
python figures/reproduce_fig3.py
python figures/reproduce_fig4_detuning.py
```

These commands write the regenerated figure files to `paper/figures/`.

## Deprecation Warning

`python/deprecated/oracle_compact.py` is retained only for provenance. It constructs the old symmetry-completed matrix using `S=S^T` and waveguide reversal. It is not the exact 6x6 scattering matrix and must not be used as the truth source. Use `python/direct_oracle.py` instead.

## Repository Status

Prepared for GitHub upload as a private reproducibility repository. The manuscript should not claim that code is publicly available until the repository URL, branch, and commit hash have been verified after upload.
