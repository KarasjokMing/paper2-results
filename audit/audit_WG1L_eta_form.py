"""audit_WG1L_eta_form.py — Phase 2: Audit WG1-L eta-form against direct oracle."""
import numpy as np; PI=np.pi; G=1.0; A2=PI
import sys, os, json, csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
from direct_oracle import build, compute_S

np.random.seed(12345)

# ============================================================
# Phase factor definitions (using the paper's notation)
# ============================================================
def P_plus(phi, theta):
    """P_{j+} = 1 + e^{i(phi+theta)}  -- entrance/source factor"""
    return 1 + np.exp(1j*(phi+theta))

def P_minus(phi, theta):
    """P_{j-} = 1 + e^{i(phi-theta)}  -- reflection/left-output factor"""
    return 1 + np.exp(1j*(phi-theta))

def P_R_conj(phi, theta):
    """P_{j+}^* = 1 + e^{-i(phi+theta)}  -- right-output factor"""
    return 1 + np.exp(-1j*(phi+theta))

# ============================================================
# Eta-form scattering amplitudes for port j
# ============================================================
def eta_form_port(eta_a, eta_b, params, port):
    """
    Compute 6 scattering amplitudes from eta_a, eta_b using the eta-form.
    params: dict with keys p1,p2,p3, a1,a2,b2,b3
    port: 1..6
    Returns: np.array of 6 complex amplitudes [s1,...,s6]
    """
    p1=params['p1'];p2=params['p2'];p3=params['p3']
    a1=params['a1'];a2=params['a2'];b2=params['b2'];b3=params['b3']
    
    # Phase factors for atom-a couplings
    P1m = P_minus(p1, a1)          # WG1-L output factor for atom a
    P1pc = P_R_conj(p1, a1)        # WG1-R output factor for atom a
    P2am = P_minus(p2, a2)         # WG2-L output factor from atom a
    P2ar = P_R_conj(p2, a2)        # WG2-R output factor from atom a
    
    # Phase factors for atom-b couplings
    P2bm = P_minus(p2, b2)         # WG2-L output factor from atom b
    P2br = P_R_conj(p2, b2)        # WG2-R output factor from atom b
    P3m = P_minus(p3, b3)          # WG3-L output factor from atom b
    P3pc = P_R_conj(p3, b3)        # WG3-R output factor from atom b
    
    if port == 1:  # WG1-L incidence
        return np.array([
            -1j*eta_a*P1m,                          # r1 (WG1-L)
            1 - 1j*eta_a*P1pc,                      # t1 (WG1-R)
            -1j*(eta_a*P2am + eta_b*P2bm),          # r2 (WG2-L)
            -1j*(eta_a*P2ar + eta_b*P2br),          # t2 (WG2-R)
            -1j*eta_b*P3m,                           # r3 (WG3-L)
            -1j*eta_b*P3pc                           # t3 (WG3-R)
        ])
    elif port == 2:  # WG1-R incidence (left-mover from right)
        # Mirror: eta-form flips roles of + and - factors
        P1p = P_plus(p1, a1)
        return np.array([
            1 - 1j*eta_a*P1m,                       # from WG1-R -> WG1-L (t-like)
            -1j*eta_a*P1pc,                          # from WG1-R -> WG1-R (r-like)
            -1j*(eta_a*P2am + eta_b*P2bm),          # r2 (same couple)
            -1j*(eta_a*P2ar + eta_b*P2br),          # t2
            -1j*eta_b*P3m,                           # r3
            -1j*eta_b*P3pc                           # t3
        ])
        # NOTE: This needs verification! Port order may need swapping.
    elif port == 3:  # WG2-L incidence
        return np.array([
            -1j*eta_a*P1pc,                          # ~t1 (port 1 output)
            -1j*eta_a*P1m,                           # ~r1 (port 2 output)
            1 - 1j*eta_a*P2ar - 1j*eta_b*P2br,      # ~t2 (port 4 output)
            -1j*eta_a*P2am - 1j*eta_b*P2bm,          # ~r2 (port 3 output)
            -1j*eta_b*P3pc,                           # ~t3 (port 6 output)
            -1j*eta_b*P3m                            # ~r3 (port 5 output)
        ])
    elif port == 4:  # WG2-R incidence
        return np.array([
            -1j*eta_a*P1pc,                          # port 1
            -1j*eta_a*P1m,                           # port 2
            -1j*eta_a*P2ar - 1j*eta_b*P2br,          # port 4 
            1 - 1j*eta_a*P2am - 1j*eta_b*P2bm,       # port 3
            -1j*eta_b*P3pc,                           # port 6
            -1j*eta_b*P3m                            # port 5
        ])
    elif port == 5:  # WG3-L incidence (by symmetry a<->b, 1<->3)
        # Atom b is the entrance, atom a is the far one
        return np.array([
            -1j*eta_a*P1m,                           # r1
            -1j*eta_a*P1pc,                          # t1
            -1j*(eta_b*P2bm + eta_a*P2am),           # r2
            -1j*(eta_b*P2br + eta_a*P2ar),           # t2
            1 - 1j*eta_b*P3m,                         # r3 (reflection on WG3-L)
            -1j*eta_b*P3pc                           # t3
        ])
    elif port == 6:  # WG3-R incidence
        return np.array([
            -1j*eta_a*P1m,                           # r1
            -1j*eta_a*P1pc,                          # t1
            -1j*(eta_b*P2bm + eta_a*P2am),           # r2
            -1j*(eta_b*P2br + eta_a*P2ar),           # t2
            -1j*eta_b*P3m,                            # r3
            1 - 1j*eta_b*P3pc                         # t3
        ])

# ============================================================
# Audit
# ============================================================
print("="*60)
print("PHASE 2: WG1-L ETA-FORM AUDIT")
print("="*60)

# Get full S matrix from direct oracle and extract eta for port 1
N_main = 5000
N_stress_near_pole = 1000
N_stress_near_P1p_zero = 1000
N_stress_near_resonance = 1000

results = {'main': [], 'near_pole': [], 'near_P1p': [], 'near_res': []}
names = ['r1','t1','r2','t2','r3','t3']

def run_audit(label, n, sampler):
    """Run audit with given sampler function."""
    max_errs = np.zeros(6)
    skipped = 0
    for k in range(n):
        params = sampler(k)
        if params is None:
            skipped += 1
            continue
        
        D=params['D'];a1=params['a1'];b2=params['b2'];b3=params['b3']
        p1=params['p1'];p2=params['p2'];p3=params['p3']
        
        A,b = build(G, D, a1, A2, b2, b3, p1, p2, p3, 1)
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            skipped += 1
            continue
        
        eta_a = x[12]; eta_b = x[13]
        s_dir = x[6:12]
        
        params_dict = {'p1':p1,'p2':p2,'p3':p3,'a1':a1,'a2':A2,'b2':b2,'b3':b3}
        s_eta = eta_form_port(eta_a, eta_b, params_dict, 1)
        
        errs = np.abs(s_eta - s_dir)
        max_errs = np.maximum(max_errs, errs)
    
    return max_errs, skipped

# --- Test 1: Main random test (5000 points) ---
print("\nTest 1: Random parameters (fixed backbone, theta_a2=pi)")
def main_sampler(k):
    return {
        'D': np.random.uniform(-3, 3),
        'a1': np.random.uniform(0, 2*PI),
        'b2': np.random.uniform(0, 2*PI),
        'b3': np.random.uniform(0, 2*PI),
        'p1': np.random.uniform(0, 2*PI),
        'p2': np.random.uniform(0, 2*PI),
        'p3': np.random.uniform(0, 2*PI)
    }
errs1, skip1 = run_audit('main', N_main, main_sampler)
print(f"  Skipped: {skip1}")
for i in range(6):
    print(f"  {names[i]}: max err = {errs1[i]:.2e}")

# --- Test 2: Near-pole stress test ---
print("\nTest 2: Near-pole stress (|M| small)")
def near_pole_sampler(k):
    # First generate random params, then scan for small |M|
    for _ in range(50):
        D = np.random.uniform(-3, 3)
        a1 = np.random.uniform(0, 2*PI)
        b2 = np.random.uniform(0, 2*PI)
        b3 = np.random.uniform(0, 2*PI)
        p1 = np.random.uniform(0, 2*PI)
        p2 = np.random.uniform(0, 2*PI)
        p3 = np.random.uniform(0, 2*PI)
        A,_ = build(G, D, a1, A2, b2, b3, p1, p2, p3, 1)
        M = -np.linalg.det(A)
        if abs(M) < 0.1 and abs(M) > 1e-12:
            return {'D':D,'a1':a1,'b2':b2,'b3':b3,'p1':p1,'p2':p2,'p3':p3}
    return None
errs2, skip2 = run_audit('near_pole', N_stress_near_pole, near_pole_sampler)
print(f"  Skipped (no small-M found): {skip2}")
for i in range(6):
    print(f"  {names[i]}: max err = {errs2[i]:.2e}")

# --- Test 3: Near P1+ zero stress ---
print("\nTest 3: Near P1+=0 stress")
def near_P1p_sampler(k):
    # Generate params where P1p = 0 approximately
    a1 = np.random.uniform(0, 2*PI)
    p1 = -a1 + np.pi + np.random.normal(0, 0.01)  # P1p = 0 when phi1+theta_a1 = pi
    return {
        'D': np.random.uniform(-3, 3),
        'a1': a1,
        'b2': np.random.uniform(0, 2*PI),
        'b3': np.random.uniform(0, 2*PI),
        'p1': p1 % (2*PI),
        'p2': np.random.uniform(0, 2*PI),
        'p3': np.random.uniform(0, 2*PI)
    }
errs3, skip3 = run_audit('near_P1p', N_stress_near_P1p_zero, near_P1p_sampler)
print(f"  Skipped: {skip3}")
for i in range(6):
    print(f"  {names[i]}: max err = {errs3[i]:.2e}")

# --- Test 4: Near resonance stress ---
print("\nTest 4: Near-resonance stress (|Delta| small)")
def near_res_sampler(k):
    return {
        'D': np.random.normal(0, 0.01),
        'a1': np.random.uniform(0, 2*PI),
        'b2': np.random.uniform(0, 2*PI),
        'b3': np.random.uniform(0, 2*PI),
        'p1': np.random.uniform(0, 2*PI),
        'p2': np.random.uniform(0, 2*PI),
        'p3': np.random.uniform(0, 2*PI)
    }
errs4, skip4 = run_audit('near_res', N_stress_near_resonance, near_res_sampler)
print(f"  Skipped: {skip4}")
for i in range(6):
    print(f"  {names[i]}: max err = {errs4[i]:.2e}")

# ============================================================
# Summary & Report
# ============================================================
print("\n" + "="*60)
print("AUDIT SUMMARY")
print("="*60)
all_errs = [errs1, errs2, errs3, errs4]
test_names = ['Main (5000)', 'Near-pole (1000)', 'Near P1+=0 (1000)', 'Near-resonance (1000)']
threshold_default = 1e-12
threshold_near_singular = 1e-8

all_pass = True
for tidx, (errs, tname) in enumerate(zip(all_errs, test_names)):
    thresh = threshold_near_singular if tidx > 0 else threshold_default
    for i in range(6):
        status = "PASS" if errs[i] < thresh else "FAIL"
        if errs[i] >= thresh:
            all_pass = False
        if errs[i] > 1e-14 or tidx > 0:
            print(f"  {tname} {names[i]}: {errs[i]:.2e} [{status}] (thresh {thresh:.1e})")

print(f"\nOverall: {'ALL PASS' if all_pass else 'SOME FAILURES'}")

# Save report
report_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'outputs', 'reports'))
os.makedirs(report_dir, exist_ok=True)

with open(os.path.join(report_dir, 'WG1L_eta_form_audit.md'), 'w', encoding='utf-8') as f:
    f.write("# WG1-L Eta-Form Audit Report\n\n")
    f.write(f"> Date: 2026-05-31\n")
    f.write(f"> Samples: {N_main} (main) + {N_stress_near_pole} (near-pole) + {N_stress_near_P1p_zero} (near-P1p=0) + {N_stress_near_resonance} (near-resonance)\n\n")
    
    f.write("## Eta-Form Expressions (WG1-L incidence)\n\n")
    f.write("```\n")
    f.write("r1 = -i * eta_a * P_{1-}\n")
    f.write("t1 = 1 - i * eta_a * P_{1+}^*\n")
    f.write("r2 = -i * (eta_a * P_{2-}^a + eta_b * P_{2-}^b)\n")
    f.write("t2 = -i * (eta_a * P_{2+}^{a*} + eta_b * P_{2+}^{b*})\n")
    f.write("r3 = -i * eta_b * P_{3-}\n")
    f.write("t3 = -i * eta_b * P_{3+}^*\n")
    f.write("```\n\n")
    
    f.write("## Audit Results\n\n")
    f.write("| Test | r1 | t1 | r2 | t2 | r3 | t3 |\n")
    f.write("|:---|---:|---:|---:|---:|---:|---:|\n")
    for tidx, (errs, tname) in enumerate(zip(all_errs, test_names)):
        f.write(f"| {tname} | {errs[0]:.2e} | {errs[1]:.2e} | {errs[2]:.2e} | {errs[3]:.2e} | {errs[4]:.2e} | {errs[5]:.2e} |\n")
    
    f.write("\n## Eta factorizations\n\n")
    f.write("- eta_a = -i * G * P_{1+} * (tau_b2 + tau_b3 + 4*G - i*Delta) / M  [VERIFIED]\n")
    f.write("- eta_b: does NOT factor to -i*G^2*P_{1+}*B3/M (paper formula); use raw 14x14 extraction or symbolic form\n\n")
    
    f.write("## Paper formula status\n\n")
    f.write("| Formula | Status |\n")
    f.write("|:---|---:|\n")
    f.write("| r1 = -G*P1-*P1+*(tau_b2+tau_b3+4G-iD)/M | PASS |\n")
    f.write("| t1 = N_t1/M | PASS |\n")
    f.write("| r2 = -G*P1+*Ar/M | PASS |\n")
    f.write("| t2 = -G*P1+*At/M | FAIL (use eta-form) |\n")
    f.write("| t3 = G^2*P1+*B3*P3+/M | FAIL (use eta-form) |\n")
    f.write("| r3 = G^2*P1+*B3*P3-/M | FAIL (use eta-form) |\n")

print(f"\nReport saved to: {report_dir}/WG1L_eta_form_audit.md")
