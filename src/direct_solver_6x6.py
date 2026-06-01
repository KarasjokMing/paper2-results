"""direct_oracle.py — True exact 6x6 S-matrix via independent 6-column solve."""
import numpy as np;PI=np.pi

def build(Gamma,Delta,thA1,thA2,thB2,thB3,phi1,phi2,phi3,port):
    A=np.zeros((14,14),dtype=complex);b=np.zeros(14,dtype=complex)
    G=Gamma;D=Delta
    e1=np.exp(1j*phi1);e1m=np.exp(-1j*phi1);e2=np.exp(1j*phi2);e2m=np.exp(-1j*phi2)
    e3=np.exp(1j*phi3);e3m=np.exp(-1j*phi3)
    sR1=1.0 if port==1 else 0.0;sL1=1.0 if port==2 else 0.0
    sR2=1.0 if port==3 else 0.0;sL2=1.0 if port==4 else 0.0
    sR3=1.0 if port==5 else 0.0;sL3=1.0 if port==6 else 0.0
    A[0,0]=-1j;A[0,12]=1.0;b[0]=-1j*sR1
    A[1,1]=1j;A[1,6]=-1j;A[1,12]=1.0
    A[2,0]=1j*e1;A[2,7]=-1j*e1;A[2,12]=np.exp(-1j*thA1)
    A[3,1]=-1j*e1m;A[3,12]=np.exp(-1j*thA1);b[3]=-1j*e1m*sL1
    A[4,2]=-1j;A[4,12]=1.0;A[4,13]=1.0;b[4]=-1j*sR2
    A[5,3]=1j;A[5,8]=-1j;A[5,12]=1.0;A[5,13]=1.0
    A[6,2]=1j*e2;A[6,9]=-1j*e2;A[6,12]=np.exp(-1j*thA2);A[6,13]=np.exp(-1j*thB2)
    A[7,3]=-1j*e2m;A[7,12]=np.exp(-1j*thA2);A[7,13]=np.exp(-1j*thB2);b[7]=-1j*e2m*sL2
    A[8,4]=-1j;A[8,13]=1.0;b[8]=-1j*sR3
    A[9,5]=1j;A[9,10]=-1j;A[9,13]=1.0
    A[10,4]=1j*e3;A[10,11]=-1j*e3;A[10,13]=np.exp(-1j*thB3)
    A[11,5]=-1j*e3m;A[11,13]=np.exp(-1j*thB3);b[11]=-1j*e3m*sL3
    A[12,0]=G*np.exp(1j*thA1)*e1;A[12,1]=G*np.exp(1j*thA1)*e1m
    A[12,2]=G*np.exp(1j*thA2)*e2;A[12,3]=G*np.exp(1j*thA2)*e2m
    A[12,6]=G;A[12,8]=G;A[12,12]=-D;b[12]=-(G*sR1+G*sR2)
    A[13,2]=G*np.exp(1j*thB2)*e2;A[13,3]=G*np.exp(1j*thB2)*e2m
    A[13,4]=G*np.exp(1j*thB3)*e3;A[13,5]=G*np.exp(1j*thB3)*e3m
    A[13,8]=G;A[13,10]=G;A[13,13]=-D;b[13]=-(G*sR2+G*sR3)
    return A,b

def compute_S(Gamma,Delta,thA1,thA2,thB2,thB3,phi1,phi2,phi3):
    """Full 6x6 S-matrix: all 6 columns independently solved."""
    S=np.zeros((6,6),dtype=complex)
    for port in range(1,7):
        A,b=build(Gamma,Delta,thA1,thA2,thB2,thB3,phi1,phi2,phi3,port)
        S[:,port-1]=np.linalg.solve(A,b)[6:12]
    return S

# Self-test
if __name__=='__main__':
    np.random.seed(42)
    max_cons=0;max_unit=0;max_TR=0
    for k in range(2000):
        D=np.random.uniform(-3,3);a1=np.random.uniform(0,2*PI)
        b2=np.random.uniform(0,2*PI);b3=np.random.uniform(0,2*PI)
        p1=np.random.uniform(0,2*PI);p2=np.random.uniform(0,2*PI);p3=np.random.uniform(0,2*PI)
        S=compute_S(1,D,a1,PI,b2,b3,p1,p2,p3)
        T=np.abs(S)**2
        mc=np.max(np.abs(T.sum(0)-1));mu=np.max(np.abs(S.conj().T@S-np.eye(6)))
        if mc>max_cons:max_cons=mc
        if mu>max_unit:max_unit=mu
        Sp=compute_S(1,D,a1,PI,b2,b3,p1,p2,p3)
        Sm=compute_S(1,D,-a1,-PI,-b2,-b3,p1,p2,p3)
        mtr=np.max(np.abs(np.abs(Sp)**2-np.abs(Sm).T**2))
        if mtr>max_TR:max_TR=mtr
    print(f"Conserv: {max_cons:.2e}  Unitarity: {max_unit:.2e}  TR: {max_TR:.2e}")
