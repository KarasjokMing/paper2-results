"""audit_all_ports_eta_form.py — Phase 4: Verify eta-form for all 6 ports."""
import numpy as np; PI=np.pi; G=1.0; A2=PI
import sys, os, json, csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
from direct_oracle import build, compute_S

def P_plus(phi, theta):
    return 1 + np.exp(1j*(phi+theta))

def P_minus(phi, theta):
    return 1 + np.exp(1j*(phi-theta))

def P_R_conj(phi, theta):
    return 1 + np.exp(-1j*(phi+theta))

# ============================================================
# Derive eta-form for each port by analyzing the jump conditions
# ============================================================
def eta_form_j(port, eta_a, eta_b, a1, a2, b2, b3, p1, p2, p3):
    """
    Eta-form expressions for port j.
    Derived from jump conditions.
    Returns 6-vector of scattering amplitudes.
    """
    # Phase factors
    P1m = P_minus(p1, a1); P1pc = P_R_conj(p1, a1)
    P1p = P_plus(p1, a1); P1mc = np.conj(P1m)  # not commonly used
    
    P2am = P_minus(p2, a2); P2ar = P_R_conj(p2, a2)
    P2bm = P_minus(p2, b2); P2br = P_R_conj(p2, b2)
    P2ap = P_plus(p2, a2); P2bp = P_plus(p2, b2)
    
    P3m = P_minus(p3, b3); P3pc = P_R_conj(p3, b3)
    P3p = P_plus(p3, b3); P3mc = np.conj(P3m)
    
    # Also need: P_{j,m}^{src} = 1 + e^{i(phi_j+theta_mj)} = P_plus(phi_j, theta_mj)
    P1ap = P_plus(p1, a1)   # = P1p
    P2ap = P_plus(p2, a2)   # = P2ap
    P2bp = P_plus(p2, b2)   # = P2bp
    P3bp = P_plus(p3, b3)   # = P3p
    
    if port == 1:  # WG1-L: right-mover incident from left
        return np.array([
            -1j*eta_a*P1m,                          # r1: WG1 left-output
            1 - 1j*eta_a*P1pc,                      # t1: WG1 right-output
            -1j*(eta_a*P2am + eta_b*P2bm),          # r2: WG2 left-output
            -1j*(eta_a*P2ar + eta_b*P2br),          # t2: WG2 right-output
            -1j*eta_b*P3m,                           # r3: WG3 left-output
            -1j*eta_b*P3pc                           # t3: WG3 right-output
        ])
    
    elif port == 2:  # WG1-R: left-mover incident from right
        # Ballistic term on WG1-L output, using left-output factor P1m
        return np.array([
            1 - 1j*eta_a*P1m,                          # WG1-L output (ballistic here)
            -1j*eta_a*P1pc,                              # WG1-R output (reflection-like)
            -1j*(eta_a*P2am + eta_b*P2bm),              # r2
            -1j*(eta_a*P2ar + eta_b*P2br),              # t2
            -1j*eta_b*P3m,                               # r3
            -1j*eta_b*P3pc                               # t3
        ])
    
    elif port == 3:  # WG2-L: right-mover incident from left of WG2
        # Paper derivation: ballistic term on WG2 right-output
        # ~t1 = -i*eta_a*P1pc, ~r1 = -i*eta_a*P1m
        # ~t2 = 1 - i*(eta_a*P2ar + eta_b*P2br)
        # ~r2 = -i*(eta_a*P2am + eta_b*P2bm)
        # ~t3 = -i*eta_b*P3pc, ~r3 = -i*eta_b*P3m
        # Port mapping: [r1,t1,r2,t2,r3,t3] -> [~r1_from_WG2, ~t1_from_WG2, ~r2, ~t2, ~r3, ~t3]
        # But careful: port 3 input -> outputs are (port1,port2,port3,port4,port5,port6)
        # Port 3 = WG2-L (output is left-mover on WG2)
        # ~r1 = WG1-L output, ~t1 = WG1-R output, ~r2 = WG2-L output, ~t2 = WG2-R output
        # This maps to: S[1,3]=~r1, S[2,3]=~t1, S[3,3]=~r2, S[4,3]=~t2, S[5,3]=~r3, S[6,3]=~t3
        return np.array([
            -1j*eta_a*P1m,                           # WG1-L output
            -1j*eta_a*P1pc,                          # WG1-R output
            -1j*(eta_a*P2am + eta_b*P2bm),           # WG2-L output
            1 - 1j*(eta_a*P2ar + eta_b*P2br),        # WG2-R output (ballistic here)
            -1j*eta_b*P3m,                           # WG3-L output
            -1j*eta_b*P3pc                           # WG3-R output
        ])
    
    elif port == 4:  # WG2-R: left-mover incident from right of WG2
        # Ballistic term on WG2-L output, using left-output factors
        return np.array([
            -1j*eta_a*P1m,                               # WG1-L
            -1j*eta_a*P1pc,                              # WG1-R
            1 - 1j*(eta_a*P2am + eta_b*P2bm),            # WG2-L (ballistic)
            -1j*(eta_a*P2ar + eta_b*P2br),               # WG2-R (reflection-like)
            -1j*eta_b*P3m,                               # WG3-L
            -1j*eta_b*P3pc                               # WG3-R
        ])
    
    elif port == 5:  # WG3-L: right-mover incident from left of WG3
        # Symmetric to port 1 with a<->b exchange
        # Atom b now is the entrance atom, atom a is the far one
        return np.array([
            -1j*eta_a*P1m,                           # WG1-L output
            -1j*eta_a*P1pc,                          # WG1-R output
            -1j*(eta_a*P2am + eta_b*P2bm),           # WG2-L output
            -1j*(eta_a*P2ar + eta_b*P2br),           # WG2-R output
            -1j*eta_b*P3m,                           # WG3-L output (reflection-like)
            1 - 1j*eta_b*P3pc                        # WG3-R output (ballistic)
        ])
    
    elif port == 6:  # WG3-R: left-mover incident from right of WG3
        # Ballistic term on WG3-L output, using left-output factor P3m
        return np.array([
            -1j*eta_a*P1m,                               # WG1-L
            -1j*eta_a*P1pc,                              # WG1-R
            -1j*(eta_a*P2am + eta_b*P2bm),              # WG2-L
            -1j*(eta_a*P2ar + eta_b*P2br),               # WG2-R
            1 - 1j*eta_b*P3m,                            # WG3-L (ballistic)
            -1j*eta_b*P3pc                               # WG3-R (reflection-like)
        ])

# ============================================================
# Audit all 6 ports
# ============================================================
np.random.seed(42)
N = 5000

port_names = {1:'WG1-L', 2:'WG1-R', 3:'WG2-L', 4:'WG2-R', 5:'WG3-L', 6:'WG3-R'}
amp_names = ['s1','s2','s3','s4','s5','s6']

max_errs = np.zeros((6, 6))  # [port-1, amp]
skipped = 0

for k in range(N):
    D = np.random.uniform(-3, 3)
    a1 = np.random.uniform(0, 2*PI)
    b2 = np.random.uniform(0, 2*PI)
    b3 = np.random.uniform(0, 2*PI)
    p1 = np.random.uniform(0, 2*PI)
    p2 = np.random.uniform(0, 2*PI)
    p3 = np.random.uniform(0, 2*PI)
    
    for port in range(1, 7):
        A, b = build(G, D, a1, A2, b2, b3, p1, p2, p3, port)
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            skipped += 1
            continue
        
        eta_a = x[12]; eta_b = x[13]
        s_dir = x[6:12]
        
        s_eta = eta_form_j(port, eta_a, eta_b, a1, A2, b2, b3, p1, p2, p3)
        
        errs = np.abs(s_eta - s_dir)
        max_errs[port-1] = np.maximum(max_errs[port-1], errs)

print("="*60)
print("PHASE 4: ALL 6 PORTS ETA-FORM AUDIT")
print("="*60)
print(f"Samples: {N}, Skipped: {skipped}")
print()

all_pass = True
for port in range(1, 7):
    pname = port_names[port]
    mx = max(max_errs[port-1])
    status = "PASS" if mx < 1e-12 else "FAIL"
    if mx >= 1e-12:
        all_pass = False
    print(f"Port {port} ({pname}): max err = {mx:.2e} [{status}]")
    if mx > 1e-14:
        for i in range(6):
            print(f"  {amp_names[i]}: {max_errs[port-1,i]:.2e}")

print(f"\nOverall: {'ALL PORTS PASS' if all_pass else 'SOME FAILURES'}")

# Save report
report_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'outputs', 'reports'))
os.makedirs(report_dir, exist_ok=True)

with open(os.path.join(report_dir, 'all_ports_eta_form.md'), 'w', encoding='utf-8') as f:
    f.write("# All 6 Ports Eta-Form Audit\n\n")
    f.write(f"> Samples: {N}, Backbone: theta_a2 = pi\n\n")
    f.write("| Port | Name | Max Error | Status |\n")
    f.write("|:---|---:|---:|:---|\n")
    for port in range(1, 7):
        mx = max(max_errs[port-1])
        status = "PASS" if mx < 1e-12 else "FAIL"
        f.write(f"| {port} | {port_names[port]} | {mx:.2e} | {status} |\n")
    f.write(f"\n## Conclusion\n\n")
    if all_pass:
        f.write("All 6 ports pass the eta-form audit.\n")
    else:
        f.write("Some ports fail. See above.\n")

print(f"\nReport saved to: {report_dir}/all_ports_eta_form.md")
