(* WG2 incidence: FULL simplification to tau-factorized form *)
(* Strategy: substitute tau, P1p, B etc. AFTER solving, then FullSimplify *)

ClearAll["Global`*"];

(* === 14 equations with src1/src3 === *)
vars = {AA,BB,CC,DD,EE,FF, t1,r1,t2,r2,t3,r3, etaA,etaB};
thA11=0; thA12=thA1; thA21=0; thA22=thA2;
thB21=0; thB22=thB2; thB31=0; thB32=thB3;

eq1 =-I(AA-src1)+etaA Exp[-I thA11];
eq2 = I(BB-r1)+etaA Exp[-I thA11];
eq3 =-I Exp[I phi1](t1-AA)+etaA Exp[-I thA12];
eq4 =-I Exp[-I phi1] BB+etaA Exp[-I thA12];
eq5 =-I(CC-src3)+etaA Exp[-I thA21]+etaB Exp[-I thB21];
eq6 = I(DD-r2)+etaA Exp[-I thA21]+etaB Exp[-I thB21];
eq7 =-I Exp[I phi2](t2-CC)+etaA Exp[-I thA22]+etaB Exp[-I thB22];
eq8 =-I Exp[-I phi2] DD+etaA Exp[-I thA22]+etaB Exp[-I thB22];
eq9 =-I EE+etaB Exp[-I thB31];
eq10= I(FF-r3)+etaB Exp[-I thB31];
eq11=-I Exp[I phi3](t3-EE)+etaB Exp[-I thB32];
eq12=-I Exp[-I phi3] FF+etaB Exp[-I thB32];

eq13 = Exp[I thA11](src1+r1) + Exp[I thA12](AA Exp[I phi1]+BB Exp[-I phi1])
     + Exp[I thA21](src3+r2) + Exp[I thA22](CC Exp[I phi2]+DD Exp[-I phi2]) 
     - I*Delta*etaA;  (* using G0=1, and note the atomic eq gives Delta*eta = sum*Gamma, 
                         with Gamma=1 and the sum already has factor 1 *)
eq14 = Exp[I thB21](src3+r2) + Exp[I thB22](CC Exp[I phi2]+DD Exp[-I phi2])
     + Exp[I thB31] r3 + Exp[I thB32](EE Exp[I phi3]+FF Exp[-I phi3]) 
     - I*Delta*etaB;

eqsAll = {eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14};

(* === Solve both incidences === *)
sol1 = Solve[Thread[(eqsAll/.{src1->1,src3->0})==0], vars][[1]];
sol3 = Solve[Thread[(eqsAll/.{src1->0,src3->1})==0], vars][[1]];

(* === Substitute tau and factorized structures === *)
tauA1 = 2 Exp[I phi1] Cos[thA1];
tauA2 = 2 Exp[I phi2] Cos[thA2];
tauB2 = 2 Exp[I phi2] Cos[thB2];
tauB3 = 2 Exp[I phi3] Cos[thB3];
P1p = 1 + Exp[I(phi1+thA1)];  P1m = 1 + Exp[I(phi1-thA1)];
P2ap = 1 + Exp[I(phi2+thA2)]; P2am = 1 + Exp[I(phi2-thA2)];
P2bp = 1 + Exp[I(phi2+thB2)]; P2bm = 1 + Exp[I(phi2-thB2)];
P3p = 1 + Exp[I(phi3+thB3)];  P3m = 1 + Exp[I(phi3-thB3)];
Bfac = 1 + Exp[I(thA2-thB2)] + Exp[I(phi2+thA2)] + Exp[I(phi2-thB2)];
d = thA2 - thB2;

(* Denominator M *)
M0 = -2 Cos[d] - 4 Exp[2 I phi2] Sin[d/2]^2
     + tauA1 tauB2 + (tauA1+tauA2) tauB3
     + (2 - I Delta)(tauA2+tauB2) + (4 - I Delta)(tauA1+tauB3)
     + 14 - 8 I Delta - Delta^2;

(* Define the simplification rules *)
tauRules = {
  Exp[I phi1] Cos[thA1] -> tauA1/2, 
  Exp[I phi2] Cos[thA2] -> tauA2/2,
  Exp[I phi2] Cos[thB2] -> tauB2/2,
  Exp[I phi3] Cos[thB3] -> tauB3/2
};

(* === AGGRESSIVELY simplify WG1 results (sanity check) === *)
Print["=== WG1: VERIFY KNOWN FORMS ==="];
Nt1check = Simplify[t1 /. sol1];
Nr1check = Simplify[r1 /. sol1];
Nt3check = Simplify[t3 /. sol1];

(* Verify r1 matches known form *)
Print["r1*M == -P1m*P1p*(tauB2+tauB3+4-I Delta): ", 
  Simplify[Nr1check*M0 + P1m*P1p*(tauB2+tauB3+4-I Delta)] == 0];

(* Verify t3 = P1p*B*P3p/M *)
Print["t3*M/(P1p*Bfac*P3p) == 1: ",
  Simplify[Nt3check*M0/(P1p*Bfac*P3p)] == 1];

(* === NOW SIMPLIFY WG2 RESULTS === *)
Print[""];
Print["=== WG2: FULLY SIMPLIFIED FORMS ==="];
Print["(All expressed in tau, P1p, Bfac, P3p, P3m where possible)"];

(* ~t1 - should equal t2 from WG1 *)
Ntt1 = FullSimplify[(t1 /. sol3) * M0];
Print["~Nt1 = ", Ntt1];
Print["  == -P1p * P2ap * (tauB3+2-I Delta) - P1p * (1-Exp[2I phi2]) * (1-Exp[-I(thB2-thA2)]) : ",
  Simplify[Ntt1 + P1p*P2ap*(tauB3+2-I Delta) + P1p*(1-Exp[2I phi2])*(1-Exp[I(thA2-thB2)])] == 0];

(* ~r1 *)
Ntr1 = FullSimplify[(r1 /. sol3) * M0];
Print["~Nr1 = ", Ntr1];

(* ~t2 - the diagonal element, has the '1' *)
Ntt2 = FullSimplify[(t2 /. sol3 - 1) * M0];  (* remove direct transmission *)
Print["~Nt2 (without direct '1') = ", Ntt2];

(* ~r2 *)
Ntr2 = FullSimplify[(r2 /. sol3) * M0];
Print["~Nr2 = ", Ntr2];

(* ~t3 *)
Ntt3 = FullSimplify[(t3 /. sol3) * M0];
Print["~Nt3 = ", Ntt3];
Print["~Nt3 / P3p = ", FullSimplify[Ntt3/P3p]];

(* ~r3 *)
Ntr3 = FullSimplify[(r3 /. sol3) * M0];
Print["~Nr3 = ", Ntr3];
Print["~Nr3 / P3m = ", FullSimplify[Ntr3/P3m]];

(* Factorize ~Nt3: check if it's proportional to P2ap * something *)
Print[""];
Print["=== FACTORIZATION CHECKS ==="];
Print["~Nt3/P3p/P2ap = ", FullSimplify[Ntt3/P3p/P2ap]];
Print["~Nr3/P3m/P2ap = ", FullSimplify[Ntr3/P3m/P2ap]];
Print["~Nt2/P1p = ", FullSimplify[Ntt2/P1p]];
Print["~Nr2/P1p = ", FullSimplify[Ntr2/P1p]];
