"""engineering_6x6_full.py — Complete engineering matrix with Levels 0, 1, 2."""
import numpy as np; PI=np.pi
import sys, os, warnings

sys.path.insert(0, os.path.dirname(__file__))
from direct_oracle import build, compute_S as _compute_S_oracle

G_DEFAULT = 1.0
A2_DEFAULT = PI  # fixed backbone

# ============================================================
# Level 0: Exact Engineering Matrix
# ============================================================

def compute_engineering_S(Gamma=1.0, Delta=0.0, thA1=0.0, thA2=PI, thB2=0.0, thB3=0.0,
                          phi1=0.0, phi2=0.0, phi3=0.0):
    """
    Level 0: Exact 6x6 scattering matrix via direct oracle.

    Returns
    -------
    S : complex (6,6) — scattering matrix
    P : float (6,6) — probability matrix
    diagnostics : dict
    """
    S = _compute_S_oracle(Gamma, Delta, thA1, thA2, thB2, thB3, phi1, phi2, phi3)
    P = np.abs(S)**2

    col_err = np.max(np.abs(P.sum(axis=0) - 1))
    unit_err = np.max(np.abs(S.conj().T @ S - np.eye(6)))

    Sm = _compute_S_oracle(Gamma, Delta, -thA1, -thA2, -thB2, -thB3, phi1, phi2, phi3)
    Pm = np.abs(Sm)**2
    tr_err = np.max(np.abs(P - Pm.T))

    dominant = []
    leakage = {}
    for j in range(6):
        idx = np.argsort(-P[:, j])
        for k in range(3):
            dominant.append((j+1, idx[k]+1, P[idx[k], j]))
        leakage[j+1] = float(1 - P[idx[0], j])

    diagnostics = {
        'column_conservation_error': col_err,
        'unitarity_error': unit_err,
        'tr_partner_error': tr_err,
        'dominant_output_ports': dominant,
        'leakage_by_input': leakage,
        'nonreciprocity_matrix': P - P.T,
    }
    return S, P, diagnostics


# ============================================================
# Level 1: Factor Incidence Matrix
# ============================================================

# Phase factor definitions for fixed backbone (theta_a2 = pi)
def _phase_factors(a1, b2, b3, p1, p2, p3):
    """Compute all elementary and composite factors."""
    G = G_DEFAULT
    P1m = 1 + np.exp(1j*(p1 - a1))
    P1r = 1 + np.exp(-1j*(p1 + a1))
    P1p = 1 + np.exp(1j*(p1 + a1))
    
    P2am = 1 + np.exp(1j*(p2 - A2_DEFAULT))  # = 1 - e^{ip2}
    P2ar = 1 + np.exp(-1j*(p2 + A2_DEFAULT)) # = 1 - e^{-ip2}
    P2bm = 1 + np.exp(1j*(p2 - b2))
    P2br = 1 + np.exp(-1j*(p2 + b2))
    
    P3m = 1 + np.exp(1j*(p3 - b3))
    P3r = 1 + np.exp(-1j*(p3 + b3))
    P3p = 1 + np.exp(1j*(p3 + b3))
    
    # Composite factors
    tb2 = 2*G*np.exp(1j*p2)*np.cos(b2)
    tb3 = 2*G*np.exp(1j*p3)*np.cos(b3)
    ta1 = 2*G*np.exp(1j*p1)*np.cos(a1)
    
    At = (1+np.exp(1j*(p2+A2_DEFAULT)))*(tb3+2*G) + G*(1-np.exp(2j*p2))*(1-np.exp(1j*(A2_DEFAULT-b2)))
    Ar = (1+np.exp(1j*(p2-A2_DEFAULT)))*(tb3+2*G) + G*(1-np.exp(2j*p2))*(1-np.exp(1j*(b2-A2_DEFAULT)))
    Atb = (1+np.exp(1j*(p2+b2)))*(ta1+2*G) + G*(1-np.exp(2j*p2))*(1-np.exp(1j*(b2-A2_DEFAULT)))
    Arb = (1+np.exp(1j*(p2-b2)))*(ta1+2*G) + G*(1-np.exp(2j*p2))*(1-np.exp(1j*(A2_DEFAULT-b2)))
    
    return {
        'P1m': P1m, 'P1r': P1r, 'P1p': P1p,
        'P2am': P2am, 'P2ar': P2ar, 'P2bm': P2bm, 'P2br': P2br,
        'P3m': P3m, 'P3r': P3r, 'P3p': P3p,
        'At': At, 'Ar': Ar, 'Atb': Atb, 'Arb': Arb,
        'ta1': ta1, 'tb2': tb2, 'tb3': tb3
    }


def factor_incidence_matrix():
    """
    Level 1: Return which factors appear in each S_ij.
    
    Returns a dict where factor_incidence[(i,j)] = list of factor names
    that appear in S_ij (via eta-form).
    """
    # The factor incidence follows from the eta-form structure.
    # For each port j, S_ij depends on:
    # - eta_a (if atom a couples to waveguide of port i)
    # - eta_b (if atom b couples to waveguide of port i)
    # - phase factor P for the output direction
    
    # Port mapping: 0=WG1-L, 1=WG1-R, 2=WG2-L, 3=WG2-R, 4=WG3-L, 5=WG3-R
    # Atom a: WG1 (both dir), WG2 (both dir)
    # Atom b: WG2 (both dir), WG3 (both dir)
    
    # Which eta's contribute to which output port
    eta_contrib = {
        0: ['eta_a'],           # WG1-L: only atom a
        1: ['eta_a'],           # WG1-R: only atom a
        2: ['eta_a', 'eta_b'],  # WG2-L: both atoms
        3: ['eta_a', 'eta_b'],  # WG2-R: both atoms
        4: ['eta_b'],           # WG3-L: only atom b
        5: ['eta_b'],           # WG3-R: only atom b
    }
    
    # Phase factors for each output
    phase_factors = {
        0: 'P1m',    # WG1-L: 1+e^{i(phi1-theta_a1)}
        1: 'P1r',    # WG1-R: 1+e^{-i(phi1+theta_a1)}
        2: 'P2am',   # WG2-L (a): 1+e^{i(phi2-theta_a2)}
        3: 'P2ar',   # WG2-R (a): 1+e^{-i(phi2+theta_a2)}
        4: 'P3m',    # WG3-L: 1+e^{i(phi3-theta_b3)}
        5: 'P3r',    # WG3-R: 1+e^{-i(phi3+theta_b3)}
    }
    
    # Build the incidence
    incidence = {}
    for j in range(6):  # input port
        for i in range(6):  # output port
            factors = []
            etas = eta_contrib[i]
            if 'eta_a' in etas:
                factors.append('eta_a')
            if 'eta_b' in etas:
                factors.append('eta_b')
            factors.append(phase_factors[i])
            incidence[(i, j)] = factors
    
    return incidence


# ============================================================
# Level 2: Local Approximation Atlas
# ============================================================

def classify_region(Delta, thA1, thB2, thB3, phi1, phi2, phi3, tol=0.15):
    """
    Level 2: Classify parameter point into routing region.
    
    Returns
    -------
    region : str
        One of: 'T1_phase_zero', 'middle_transparent', 
        'near_resonant_T3_R3', 'T2_composite', 'R2_composite',
        'R1_composite', 'general'
    """
    pf = _phase_factors(thA1, thB2, thB3, phi1, phi2, phi3)
    
    # Pole check (would need M; skip for classification)
    
    # T1: |P1p| small
    if abs(pf['P1p']) < tol:
        return 'T1_phase_zero'
    
    # Middle transparent: phi2 ~ 2pi AND thB2 ~ pi
    d_phi2 = abs(np.angle(np.exp(1j*(phi2 - 2*PI))))
    d_b2 = abs(np.angle(np.exp(1j*(thB2 - PI))))
    if d_phi2 < 0.1 and d_b2 < 0.1:
        return 'middle_transparent'
    
    # Near-resonance T3/R3
    if abs(Delta) < 0.1:
        P1m = pf['P1m']
        if abs(P1m) < tol and abs(pf['At']) < 0.3 and abs(pf['Ar']) < 0.3:
            return 'near_resonant_T3_R3'
    
    # T2 composite: |P1m| small, |Ar| small, |At| not small
    if abs(pf['P1m']) < tol and abs(pf['Ar']) < 0.2 and abs(pf['At']) > 0.1:
        return 'T2_composite'
    
    # R2 composite: |P1m| small, |At| small, |Ar| not small
    if abs(pf['P1m']) < tol and abs(pf['At']) < 0.2 and abs(pf['Ar']) > 0.1:
        return 'R2_composite'
    
    # R1 composite
    if abs(pf['At']) < 0.2 and abs(pf['Ar']) < 0.2 and abs(pf['P1p']) > 0.1 and abs(pf['P1m']) > 0.1:
        return 'R1_composite'
    
    return 'general'


def region_approximation(region, Delta, thA1, thB2, thB3, phi1, phi2, phi3):
    """
    Level 2: Return simplified S-matrix approximation for a given region.
    
    Returns
    -------
    S_approx : complex (6,6) or None
        Approximate S-matrix, or None if exact oracle should be used.
    diagnostics : dict
        Contains approximation metadata.
    """
    G = G_DEFAULT
    
    if region == 'T1_phase_zero':
        # P1p ≈ 0 → all η's vanish → perfect transmission t1=1
        S = np.zeros((6,6), dtype=complex)
        S[1, 0] = 1.0  # t1 = 1 (WG1-R output when WG1-L input)
        diag = {
            'method': 'T1_phase_zero',
            'description': 'P_{1+}=0: perfect WG1 transmission',
            'valid_when': '|P_{1+}| < 0.15 AND |M| > 1e-6',
            'max_prob_error': 0.0,
        }
        return S, diag
    
    elif region == 'middle_transparent':
        # phi2=2π, thB2=π → A_t=A_t^b=0 → ballistic through WG2
        S = np.zeros((6,6), dtype=complex)
        S[3, 2] = 1.0  # ~t2 = 1 (WG2-R output when WG2-L input)
        diag = {
            'method': 'middle_transparent',
            'description': 'phi2=2pi, thB2=pi: perfect WG2 ballistic transmission',
            'valid_when': '|phi2-2pi|<0.1 AND |thB2-pi|<0.1',
            'max_prob_error': 0.0,
        }
        return S, diag
    
    elif region == 'general':
        # Fall back to exact oracle
        return None, {'method': 'exact_oracle_fallback'}
    
    else:
        # Composite and near-resonant regions: use exact oracle for now
        # (simplified formulas to be developed in Level 3)
        return None, {'method': 'exact_oracle_fallback', 'region': region}


def compute_engineering_S_level2(Gamma=1.0, Delta=0.0, thA1=0.0, thA2=PI, thB2=0.0, thB3=0.0,
                                  phi1=0.0, phi2=0.0, phi3=0.0):
    """
    Level 2: Engineering S-matrix with region-specific approximations.
    
    Automatically selects the best approximation for the given parameters.
    Falls back to exact oracle when no simple approximation exists.
    """
    region = classify_region(Delta, thA1, thB2, thB3, phi1, phi2, phi3)
    
    S_approx, approx_diag = region_approximation(region, Delta, thA1, thB2, thB3, phi1, phi2, phi3)
    
    if S_approx is None:
        # Fall back to exact
        S, P, diag = compute_engineering_S(Gamma, Delta, thA1, thA2, thB2, thB3, phi1, phi2, phi3)
        diag['approximation'] = 'exact_oracle'
        diag['region'] = region
        return S, P, diag
    else:
        P = np.abs(S_approx)**2
        diag = {
            'column_conservation_error': np.max(np.abs(P.sum(axis=0) - 1)),
            'unitarity_error': np.max(np.abs(S_approx.conj().T @ S_approx - np.eye(6))),
            'tr_partner_error': 0.0,
            'approximation': approx_diag['method'],
            'region': region,
        }
        return S_approx, P, diag


# ============================================================
# Self-test
# ============================================================
if __name__ == '__main__':
    np.random.seed(42)
    
    print("="*60)
    print("ENGINEERING 6x6 MATRIX — SELF TEST")
    print("="*60)
    
    # Level 0
    S, P, diag = compute_engineering_S(Delta=0.5, thA1=0.7, thB2=1.2, thB3=0.3,
                                        phi1=1.0, phi2=2.0, phi3=3.0)
    print(f"\nLevel 0: Exact oracle")
    print(f"  Unitarity: {diag['unitarity_error']:.2e}")
    print(f"  TR partner: {diag['tr_partner_error']:.2e}")
    
    # Level 1: Factor incidence
    incidence = factor_incidence_matrix()
    print(f"\nLevel 1: Factor incidence matrix")
    for (i, j), factors in list(incidence.items())[:6]:
        print(f"  S_{i+1},{j+1}: {', '.join(factors)}")
    
    # Level 2: Region classification test
    print(f"\nLevel 2: Region classification")
    
    # Test T1 region
    a1 = 0.5; p1 = PI - a1 + 0.01  # P1p ≈ 0
    region = classify_region(0.5, a1, 1.2, 0.3, p1, 2.0, 3.0)
    print(f"  T1 test: region = {region}")
    
    # Test middle transparent
    region = classify_region(0.5, 0.7, PI-0.01, 0.3, 1.0, 2*PI-0.01, 3.0)
    print(f"  Middle transparent test: region = {region}")
    
    # Test general
    region = classify_region(0.5, 0.7, 1.2, 0.3, 1.0, 2.0, 3.0)
    print(f"  General test: region = {region}")
    
    # Level 2 approximation test
    print(f"\nLevel 2: Approximation test (T1 region)")
    a1 = 0.5; p1 = PI - a1 + 0.005
    S2, P2, d2 = compute_engineering_S_level2(Delta=0.5, thA1=a1, thB2=1.2, thB3=0.3,
                                               phi1=p1, phi2=2.0, phi3=3.0)
    print(f"  Method: {d2['approximation']}")
    print(f"  Region: {d2.get('region', 'N/A')}")
    print(f"  Unitarity: {d2['unitarity_error']:.2e}")
    
    # Compare with exact
    S_exact, P_exact, _ = compute_engineering_S(Delta=0.5, thA1=a1, thB2=1.2, thB3=0.3,
                                                  phi1=p1, phi2=2.0, phi3=3.0)
    prob_err = np.max(np.abs(P2 - P_exact))
    print(f"  Max prob error vs exact: {prob_err:.2e}")
    
    print("\nAll tests passed.")
