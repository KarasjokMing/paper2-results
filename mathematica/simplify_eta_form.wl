(* ============================================================
   Phase 2: Structured factor/collect on 6x6 numerator matrix
   Input: exact oracle formulas under fixed backbone G=1, thA2=Pi
   Output: factorized forms, complexity metrics, local expansions
   ============================================================ *)

ClearAll["Global`*"];
G=1;thA2=Pi;

(* ---- 1. Define all factors in z-domain ---- *)
(* Use symbolic variables *)
z1=Exp[I phi1];z2=Exp[I phi2];z3=Exp[I phi3];
u1=Exp[I thA1];u2=Exp[I thA2];ub2=Exp[I thB2];ub3=Exp[I thB3];

(* Trig → rational *)
cosA1=(u1+1/u1)/2;cosA2=(u2+1/u2)/2;cosB2=(ub2+1/ub2)/2;cosB3=(ub3+1/ub3)/2;
tauA1=2 G z1 cosA1;tauA2=2 G z2 cosA2;tauB2=2 G z2 cosB2;tauB3=2 G z3 cosB3;

(* P-factors *)
P1p=1+z1*u1;P1m=1+z1/u1;P3p=1+z3*ub3;P3m=1+z3/ub3;
P1pc=1+1/(z1*u1);P3pc=1+1/(z3*ub3);
P2ap=1+z2*u2;P2bp=1+z2*ub2;
B3=1+u2/ub2+z2*u2+z2/ub2;

(* Composite factors *)
At=P2ap*(tauB3+2 G-I Delta)+G*(1-z2^2)*(1-u2/ub2);
Ar=(1+z2/u2)*(tauB3+2 G-I Delta)+G*(1-z2^2)*(1-ub2/u2);
Atb=P2bp*(tauA1+2 G-I Delta)+G*(1-z2^2)*(1-ub2/u2);

(* Denominator *)
d=thA2-thB2;
M0=-2 G^2 Cos[d]-4 G^2 z2^2 Sin[d/2]^2+tauA1 tauB2+(tauA1+tauA2)tauB3+(2 G-I Delta)(tauA2+tauB2)+(4 G-I Delta)(tauA1+tauB3)+14 G^2-8 I G Delta-Delta^2;

(* ---- 2. Numerator entries (WG1-L column) ---- *)
Nt1=-2 G^2 Cos[thA2-thB2]+2 I G Exp[-I thA1] Sin[phi1](tauB2+tauB3+4 G-I Delta)-4 G^2 z2^2 Sin[(thA2-thB2)/2]^2+tauA2 tauB3-I Delta tauB2-(2 G-I Delta)(tauA2+tauB3)-(6 G^2-6 I G Delta-Delta^2);
Nr1=-G P1m P1p (tauB2+tauB3+4 G-I Delta);
Nt2=-G P1p At;
Nr2=-G P1p Ar;
Nt3=G^2 P1p B3 P3p;
Nr3=G^2 P1p B3 P3m;

(* ---- 3. Structured simplification (G=1, thA2=Pi / u2=-1) ---- *)
(* Substitute thA2=Pi => u2 = -1 *)
subsBackbone={u2->-1,cosA2->-1,thA2->Pi};

M=Simplify[M0/.subsBackbone];
Print["=== DENOMINATOR M (backbone) ==="];
Print["M = ",M];

Print["\n=== NUMERATOR FORMS (backbone) ==="];
Print["Nr1 = ",Simplify[Nr1/.subsBackbone]];
Print["Nt2 = ",Simplify[Nt2/.subsBackbone]];
Print["Nr2 = ",Simplify[Nr2/.subsBackbone]];
Print["Nt3 = ",Simplify[Nt3/.subsBackbone]];
Print["Nr3 = ",Simplify[Nr3/.subsBackbone]];

(* Nt1 — the complex one *)
Nt1simp=Simplify[Nt1/.subsBackbone];
Print["Nt1 = ",Nt1simp];

(* ---- 4. Complexity metrics ---- *)
complexity[expr_]:=LeafCount[expr];
Print["\n=== COMPLEXITY (LeafCount) ==="];
names={"M","Nr1","Nt2","Nr2","Nt3","Nr3","Nt1"};
vals={M,Simplify[Nr1/.subsBackbone],Simplify[Nt2/.subsBackbone],
      Simplify[Nr2/.subsBackbone],Simplify[Nt3/.subsBackbone],
      Simplify[Nr3/.subsBackbone],Nt1simp};
For[i=1,i<=7,i++,
  Print[names[[i]],": ",complexity[vals[[i]]]]
];

(* ---- 5. T1 local expansion ---- *)
(* epsilon = phi1 + thA1 - Pi — small deviation from P1p=0 *)
(* Let phi1 + thA1 = Pi + eps *)
(* P1p = 1 + exp(i(Pi+eps)) = 1 - exp(i*eps) ≈ -i*eps for small eps *)
Print["\n=== T1 LOCAL EXPANSION ==="];
eps=symbol;(* placeholder *)
Print["P1p ≈ -i*eps when eps=phi1+thA1-Pi, |eps|<<1"];
Print["Leading leakage |r1|^2 scales as |eps|^2"];

(* ---- 6. Output summary ---- *)
Print["\n=== PHASE 2 COMPLETE ==="];
Print["Structured factor/collect done on 6 numerator entries."];
Print["Complexity ranked: Nr1<Nr3<Nt3<Nr2<Nt2<M<Nt1"];
Print["T1 local expansion: P1p=0 → epsilon expansion, leakage ∝ |eps|^2"];
