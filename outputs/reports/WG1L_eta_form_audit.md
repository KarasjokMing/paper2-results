# WG1-L Eta-Form Audit Report

> Date: 2026-05-31
> Samples: 5000 (main) + 1000 (near-pole) + 1000 (near-P1p=0) + 1000 (near-resonance)

## Eta-Form Expressions (WG1-L incidence)

```
r1 = -i * eta_a * P_{1-}
t1 = 1 - i * eta_a * P_{1+}^*
r2 = -i * (eta_a * P_{2-}^a + eta_b * P_{2-}^b)
t2 = -i * (eta_a * P_{2+}^{a*} + eta_b * P_{2+}^{b*})
r3 = -i * eta_b * P_{3-}
t3 = -i * eta_b * P_{3+}^*
```

## Audit Results

| Test | r1 | t1 | r2 | t2 | r3 | t3 |
|:---|---:|---:|---:|---:|---:|---:|
| Main (5000) | 1.23e-15 | 1.30e-15 | 1.18e-15 | 1.67e-15 | 1.01e-15 | 1.09e-15 |
| Near-pole (1000) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| Near P1+=0 (1000) | 7.23e-16 | 5.58e-16 | 4.48e-16 | 3.78e-16 | 3.58e-16 | 3.33e-16 |
| Near-resonance (1000) | 1.07e-15 | 9.49e-16 | 1.22e-15 | 1.56e-15 | 9.61e-16 | 8.78e-16 |

## Eta factorizations

- eta_a = -i * G * P_{1+} * (tau_b2 + tau_b3 + 4*G - i*Delta) / M  [VERIFIED]
- eta_b: does NOT factor to -i*G^2*P_{1+}*B3/M (paper formula); use raw 14x14 extraction or symbolic form

## Paper formula status

| Formula | Status |
|:---|---:|
| r1 = -G*P1-*P1+*(tau_b2+tau_b3+4G-iD)/M | PASS |
| t1 = N_t1/M | PASS |
| r2 = -G*P1+*Ar/M | PASS |
| t2 = -G*P1+*At/M | FAIL (use eta-form) |
| t3 = G^2*P1+*B3*P3+/M | FAIL (use eta-form) |
| r3 = G^2*P1+*B3*P3-/M | FAIL (use eta-form) |
