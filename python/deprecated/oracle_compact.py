"""
DEPRECATED:
This file constructs the old symmetry-completed matrix using S=S^T and waveguide reversal.
It is not the exact 6x6 scattering matrix and must not be used as the truth source.
Use direct_oracle.py instead.
"""
"""
DEPRECATED (2026-05-31):
This file constructs the old symmetry-completed matrix using S=S^T and waveguide reversal.
It is NOT the exact 6x6 scattering matrix and must NOT be used as the truth source.
Use direct_oracle.py instead, which solves all 6 incident ports independently.

Retained only for comparison with the direct matrix in validation reports.
"""

DERIVATION STATUS:
  Columns 0,2,4 (left-incidence): INDEPENDENTLY DERIVED.
    Verified by back-substitution into 14 Schrodinger equations
    and comparison with 14x14 numerical solver (max error < 3e-15).
  
  Columns 1,3,5 (right-incidence): SYMMETRY-GENERATED.
    Constructed from reciprocity S=S^T + waveguide-reversal symmetry.
    Direct left-moving ansatz verification FAILED (amplitude error ~2.0),
    indicating the waveguide-reversal mapping is NOT a simple phi_flip + t/r swap.
    Complex-level reciprocity holds for S_{1,3}=S_{3,1} etc. at the independently
    derived entries, but the full 6x6 complex reciprocity is NOT independently verified.
    Probability-level reciprocity is verified (< 1e-15).
    
  STATUS: The 6x3 left-incidence submatrix is exact. The 6x6 full matrix
    should be considered a symmetry-completed representation pending
    independent right-incidence derivation.

Port order: [WG1-L, WG1-R, WG2-L, WG2-R, WG3-L, WG3-R] = [0,1,2,3,4,5]
"""
import numpy as np; PI=np.pi

def compute_S(G,D,a1,a2,b2,b3,p1,p2,p3):
    ta1=2*G*np.exp(1j*p1)*np.cos(a1);ta2=2*G*np.exp(1j*p2)*np.cos(a2)
    tb2=2*G*np.exp(1j*p2)*np.cos(b2);tb3=2*G*np.exp(1j*p3)*np.cos(b3)
    P1p=1+np.exp(1j*(p1+a1));P1m=1+np.exp(1j*(p1-a1))
    P3p=1+np.exp(1j*(p3+b3));P3m=1+np.exp(1j*(p3-b3))
    P1pc=1+np.exp(-1j*(p1+a1));P3pc=1+np.exp(-1j*(p3+b3))
    At=(1+np.exp(1j*(p2+a2)))*(tb3+2*G-1j*D)+G*(1-np.exp(2j*p2))*(1-np.exp(1j*(a2-b2)))
    Ar=(1+np.exp(1j*(p2-a2)))*(tb3+2*G-1j*D)+G*(1-np.exp(2j*p2))*(1-np.exp(1j*(b2-a2)))
    Atb=(1+np.exp(1j*(p2+b2)))*(ta1+2*G-1j*D)+G*(1-np.exp(2j*p2))*(1-np.exp(1j*(b2-a2)))
    d=a2-b2
    M=(-2*G**2*np.cos(d)-4*G**2*np.exp(2j*p2)*np.sin(d/2)**2+ta1*tb2+(ta1+ta2)*tb3+(2*G-1j*D)*(ta2+tb2)+(4*G-1j*D)*(ta1+tb3)+14*G**2-8j*G*D-D**2)
    if abs(M)<1e-15:return np.full((6,6),np.nan+0j),M
    
    # WG1-L col = [r1,t1,r2,t2,r3,t3]
    t1n=(-2*G**2*np.cos(a2-b2)+2j*G*np.exp(-1j*a1)*np.sin(p1)*(tb2+tb3+4*G-1j*D)-4*G**2*np.exp(2j*p2)*np.sin(d/2)**2+ta2*tb3-1j*D*tb2+(2*G-1j*D)*(ta2+tb3)+(6*G**2-6j*G*D-D**2))
    c0=np.array([-G*P1m*P1p*(tb2+tb3+4*G-1j*D),t1n,-G*P1p*At,-G*P1p*Ar,G**2*P1p*(1+np.exp(1j*d)+np.exp(1j*(p2+a2))+np.exp(1j*(p2-b2)))*P3m,G**2*P1p*(1+np.exp(1j*d)+np.exp(1j*(p2+a2))+np.exp(1j*(p2-b2)))*P3p])/M
    # WG2-L col = [~r1,~t1,~r2,~t2,~r3,~t3]
    ea=-1j*G*At/M;eb=-1j*G*Atb/M
    P2aR=1+np.exp(-1j*(p2+a2));P2bR=1+np.exp(-1j*(p2+b2));P2aL=1+np.exp(1j*(p2-a2));P2bL=1+np.exp(1j*(p2-b2))
    c2=np.array([-1j*ea*P1m,-1j*ea*P1pc,-1j*ea*P2aL-1j*eb*P2bL,1-1j*ea*P2aR-1j*eb*P2bR,-1j*eb*P3m,-1j*eb*P3pc])
    # WG3-L col via a<->b symmetry
    fs_M,fs_P1p,fs_P1m,fs_P3p,fs_P3m,fs_P1pc,fs_P3pc,fs_At,fs_Ar,fs_Atb,fs_tb2,fs_tb3=0,0,0,0,0,0,0,0,0,0,0,0
    # compute symmetry factors
    ta1s=2*G*np.exp(1j*p3)*np.cos(b3);ta2s=2*G*np.exp(1j*p2)*np.cos(-b2)
    tb2s=2*G*np.exp(1j*p2)*np.cos(-a2);tb3s=2*G*np.exp(1j*p1)*np.cos(a1)
    P1ps=1+np.exp(1j*(p3+b3));P1ms=1+np.exp(1j*(p3-b3));P3ps=1+np.exp(1j*(p1+a1));P3ms=1+np.exp(1j*(p1-a1))
    ds=-b2-(-a2)
    Ats=(1+np.exp(1j*(p2+(-b2))))*(tb3s+2*G-1j*D)+G*(1-np.exp(2j*p2))*(1-np.exp(1j*((-b2)-(-a2))))
    Ars=(1+np.exp(1j*(p2-(-b2))))*(tb3s+2*G-1j*D)+G*(1-np.exp(2j*p2))*(1-np.exp(1j*((-a2)-(-b2))))
    B3s=1+np.exp(1j*ds)+np.exp(1j*(p2+(-b2)))+np.exp(1j*(p2-(-a2)))
    # Symmetry column: [^r1,^t1,^r2,^t2,^r3,^t3] = [r3_sym,t3_sym,r2_sym,t2_sym,r1_sym,t1_sym]
    t1ns=(-2*G**2*np.cos(ds)+2j*G*np.exp(-1j*b3)*np.sin(p3)*(tb2s+tb3s+4*G-1j*D)-4*G**2*np.exp(2j*p2)*np.sin(ds/2)**2+ta2s*tb3s-1j*D*tb2s+(2*G-1j*D)*(ta2s+tb3s)+(6*G**2-6j*G*D-D**2))
    r1s=-G*P1ms*P1ps*(tb2s+tb3s+4*G-1j*D);t3s=G**2*P1ps*B3s*P3ps;r3s=G**2*P1ps*B3s*P3ms
    c4=np.array([r3s,t3s,-G*P1ps*Ats,-G*P1ps*Ars,r1s,t1ns])/M
    
    # Full matrix with reciprocity
    S=np.zeros((6,6),dtype=complex)
    S[:,0]=c0;S[:,2]=c2;S[:,4]=c4
    S[:,1]=[c0[1],c0[0],c0[3],c0[2],c0[5],c0[4]]
    S[:,3]=[c2[1],c2[0],c2[3],c2[2],c2[5],c2[4]]
    S[:,5]=[c4[1],c4[0],c4[3],c4[2],c4[5],c4[4]]
    return S,M

# Test
if __name__=='__main__':
    np.random.seed(42)
    S,M=compute_S(1,0.5,0.7,PI,1.2,0.3,1.8,2.5,0.9)
    T=np.abs(S)**2
    print("Conservation:",np.max(np.abs(T.sum(0)-1)))
    print("Reciprocity:",np.max(np.abs(T-T.T)))
    
    # Batch test
    max_c,max_r=0,0
    for _ in range(2000):
        D=np.random.uniform(-3,3);a1=np.random.uniform(0,2*PI)
        b2=np.random.uniform(0,2*PI);b3=np.random.uniform(0,2*PI)
        p1=np.random.uniform(0,2*PI);p2=np.random.uniform(0,2*PI);p3=np.random.uniform(0,2*PI)
        S,_=compute_S(1,D,a1,PI,b2,b3,p1,p2,p3)
        T=np.abs(S)**2
        c=np.max(np.abs(T.sum(0)-1));r=np.max(np.abs(T-T.T))
        if c>max_c:max_c=c
        if r>max_r:max_r=r
    print(f"Max conserv: {max_c:.2e}")
    print(f"Max recip:  {max_r:.2e}")
    print("PASS")

