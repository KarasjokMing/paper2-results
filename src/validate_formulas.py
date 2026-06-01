"""validate_direct_oracle.py — Comprehensive Phase 1 validation of direct oracle."""
import numpy as np; PI=np.pi; G=1.0; A2=PI
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from direct_solver_6x6 import build, compute_S

np.random.seed(12345)
N = 5000

results = {
    'column_conservation': [],
    'unitarity': [],
    'tr_partner': [],
    'zero_phase_reciprocity': [],
    'det_M': [],
    'eta_check': []
}

# === Test 1: Fixed backbone (theta_a2 = pi) ===
print("Test 1: Fixed backbone (theta_a2 = pi), 5000 random points...")
for k in range(N):
    D = np.random.uniform(-3, 3)
    a1 = np.random.uniform(0, 2*PI)
    b2 = np.random.uniform(0, 2*PI)
    b3 = np.random.uniform(0, 2*PI)
    p1 = np.random.uniform(0, 2*PI)
    p2 = np.random.uniform(0, 2*PI)
    p3 = np.random.uniform(0, 2*PI)
    
    Sp = compute_S(G, D, a1, A2, b2, b3, p1, p2, p3)
    Tp = np.abs(Sp)**2
    
    # Column conservation
    cc = np.max(np.abs(Tp.sum(axis=0) - 1))
    results['column_conservation'].append(float(cc))
    
    # Unitarity
    uu = np.max(np.abs(Sp.conj().T @ Sp - np.eye(6)))
    results['unitarity'].append(float(uu))
    
    # TR partner
    Sm = compute_S(G, D, -a1, -A2, -b2, -b3, p1, p2, p3)
    Tm = np.abs(Sm)**2
    tr = np.max(np.abs(Tp - Tm.T))
    results['tr_partner'].append(float(tr))
    
    # M = -det(A) check
    A1, b1 = build(G, D, a1, A2, b2, b3, p1, p2, p3, 1)
    detA = np.linalg.det(A1)
    results['det_M'].append(float(abs(detA)))

# === Test 2: Zero-phase reciprocity ===
print("Test 2: Zero-phase reciprocity, 1000 random points...")
for k in range(1000):
    D = np.random.uniform(-3, 3)
    p1 = np.random.uniform(0, 2*PI)
    p2 = np.random.uniform(0, 2*PI)
    p3 = np.random.uniform(0, 2*PI)
    
    Sp = compute_S(G, D, 0, 0, 0, 0, p1, p2, p3)
    Tp = np.abs(Sp)**2
    zr = np.max(np.abs(Tp - Tp.T))
    results['zero_phase_reciprocity'].append(float(zr))

# === Summary ===
print("\n" + "="*60)
print("PHASE 1: DIRECT ORACLE VALIDATION RESULTS")
print("="*60)

thresholds = {
    'column_conservation': ('Column conservation', 1e-12),
    'unitarity': ('Unitarity ||Sdag S-I||_inf', 1e-10),
    'tr_partner': ('TR partner relation', 1e-10),
    'zero_phase_reciprocity': ('Zero-phase reciprocity', 1e-10),
}

all_pass = True
for key, (name, thresh) in thresholds.items():
    vals = results[key]
    mx = max(vals)
    passed = mx < thresh
    status = "PASS" if passed else "FAIL"
    if not passed:
        all_pass = False
    print(f"  {name:<35s}: max={mx:.2e}  [{status}] (threshold {thresh:.1e})")

print(f"\n  det(A) mean abs:         {np.mean(results['det_M']):.6e}")
print(f"  Overall: {'ALL PASS' if all_pass else 'SOME FAILURES'}")

# === Save report ===
report_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'reports'))
os.makedirs(report_dir, exist_ok=True)

with open(os.path.join(report_dir, 'validation_report.md'), 'w') as f:
    f.write("# Direct Oracle Validation Report\n\n")
    f.write(f"> Date: 2026-05-31\n")
    f.write(f"> Samples: {N} (Test 1) + 1000 (Test 2)\n\n")
    f.write("## Results\n\n")
    f.write("| Test | Max Error | Threshold | Status |\n")
    f.write("|:---|---:|---:|:---|\n")
    for key, (name, thresh) in thresholds.items():
        mx = max(results[key])
        status = "PASS" if mx < thresh else "FAIL"
        f.write(f"| {name} | {mx:.2e} | {thresh:.1e} | {status} |\n")
    f.write(f"\n## Conclusion\n\n")
    if all_pass:
        f.write("All validation tests passed. The direct oracle is a valid truth source.\n")
    else:
        f.write("Some tests failed. See above for details.\n")

print(f"\nReport saved to: {report_dir}/validation_report.md")
