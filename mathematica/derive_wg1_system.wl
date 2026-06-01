(*
  Direct 14-equation solver for WG1 and WG2 incidence
  Following the methodology of Multilevel_quantum_routing_2atom.pdf:
  write ALL jump conditions + atomic equations, Solve[] directly.
  No Cramer's rule, no intermediate eta variables.
*)

(* === PHYSICAL PARAMETERS === *)
(* Use Gamma, Delta, and phase angles *)
(* tau_mj = 2 Gamma Exp[I phi_j] Cos[theta_mj] *)

(* === UNKNOWNS === *)
(* Waveguide intermediate fields: A,B (WG1), C,D (WG2), E,F (WG3) *)
(* Scattering outputs: t1,r1,t2,r2,t3,r3 *)
(* Atomic amplitudes: etaA, etaB   (eta = g*f/v_g, dimensionless) *)

(* === EQUATIONS FOR WG1 LEFT INCIDENCE === *)
(* 
  Jump conditions use the ansatz:
  WG1: Phi_R1 = e^{ik1 x}[Theta(-x) + A Theta(x)Theta(l1-x) + t1 Theta(x-l1)]
        Phi_L1 = e^{-ik1 x}[r1 Theta(-x) + B Theta(x)Theta(l1-x)]
  WG2: Phi_R2 = e^{ik2 x}[C Theta(x)Theta(l2-x) + t2 Theta(x-l2)]
        Phi_L2 = e^{-ik2 x}[r2 Theta(-x) + D Theta(x)Theta(l2-x)]
  WG3: Phi_R3 = e^{ik3 x}[E Theta(x)Theta(l3-x) + t3 Theta(x-l3)]
        Phi_L3 = e^{-ik3 x}[r3 Theta(-x) + F Theta(x)Theta(l3-x)]

  Note: For WG1 incidence, the "1" in Phi_R1 at x=0- is the source.
  For WG2 incidence, the "1" moves to Phi_R2 at x=0-.
*)

(* Let's build the equations step by step using the 
   Appendix A2 equations from the reference paper *)

ClearAll[Gamma, Delta, phi1, phi2, phi3, thA1, thA2, thB2, thB3];
ClearAll[A, B, C, D, E, F, t1, r1, t2, r2, t3, r3, etaA, etaB];
ClearAll[thA11, thA12, thA21, thA22, thB21, thB22, thB31, thB32];

(* Individual coupling phases *)
thA11 = 0; thA12 = thA1;
thA21 = 0; thA22 = thA2;
thB21 = 0; thB22 = thB2;
thB31 = 0; thB32 = thB3;

(* Build equations for general incidence: source1 at WG1 x=0, source3 at WG2 x=0 *)
(* source1=1 for WG1 incidence, source1=0 for WG2 incidence *)
(* source3=1 for WG2 incidence, source3=0 for WG1 incidence *)

(* === WG1 COUPLING POINTS === *)
(* x=0, R-field *)
eq1 = -I (A - src1) + etaA Exp[-I thA11];
(* x=0, L-field *)
eq2 = I (B - r1) + etaA Exp[-I thA11];
(* x=l1, R-field *)
eq3 = -I Exp[I phi1] (t1 - A) + etaA Exp[-I thA12];
(* x=l1, L-field *)
eq4 = -I Exp[-I phi1] B + etaA Exp[-I thA12];

(* === WG2 COUPLING POINTS === *)
(* x=0, R-field -- source3 enters here *)
eq5 = -I (C - src3) + etaA Exp[-I thA21] + etaB Exp[-I thB21];
(* x=0, L-field *)
eq6 = I (D - r2) + etaA Exp[-I thA21] + etaB Exp[-I thB21];
(* x=l2, R-field *)
eq7 = -I Exp[I phi2] (t2 - C) + etaA Exp[-I thA22] + etaB Exp[-I thB22];
(* x=l2, L-field *)
eq8 = -I Exp[-I phi2] D + etaA Exp[-I thA22] + etaB Exp[-I thB22];

(* === WG3 COUPLING POINTS === *)
(* x=0, R-field *)
eq9 = -I E + etaB Exp[-I thB31];
(* x=0, L-field *)
eq10 = I (F - r3) + etaB Exp[-I thB31];
(* x=l3, R-field *)
eq11 = -I Exp[I phi3] (t3 - E) + etaB Exp[-I thB32];
(* x=l3, L-field *)
eq12 = -I Exp[-I phi3] F + etaB Exp[-I thB32];

(* === ATOMIC EQUATIONS === *)
(* Atom a: sum over WG1, WG2 coupling points *)
(* From Eq. (A1g) in the reference *)
(* For WG1 incidence (src1=1):
   g [e^{i thA11}(Phi_R1(0)+Phi_L1(0)) + e^{i thA12}(Phi_R1(l1)+Phi_L1(l1))
      + e^{i thA21}(Phi_R2(0)+Phi_L2(0)) + e^{i thA22}(Phi_R2(l2)+Phi_L2(l2))] 
   = Delta * f_a^g
*)
(* In eta variables: divide by v_g, use eta = g f / v_g, Gamma = g^2/v_g *)
(* The atomic equations from the reference are equivalent to:
   A11 etaA + A12 etaB = S_a
   A21 etaA + A22 etaB = S_b
   But let's build them directly from the field ansatz. *)

(* Fields at coupling points (using the ansatz): *)
(* WG1 x=0: Phi_R1(0)=src1, Phi_L1(0)=r1 *)
(* WG1 x=l1: Phi_R1(l1)=A Exp[I phi1], Phi_L1(l1)=B Exp[-I phi1] *)
(* WG2 x=0: Phi_R2(0)=src3, Phi_L2(0)=r2 *)
(* WG2 x=l2: Phi_R2(l2)=C Exp[I phi2], Phi_L2(l2)=D Exp[-I phi2] *)
(* WG3 x=0: Phi_R3(0)=0, Phi_L3(0)=r3 *)
(* WG3 x=l3: Phi_R3(l3)=E Exp[I phi3], Phi_L3(l3)=F Exp[-I phi3] *)

(* Atom a equation: *)
(* g [e^{i thA11}(src1+r1) + e^{i thA12}(A e^{i phi1} + B e^{-i phi1})
      + e^{i thA21}(src3+r2) + e^{i thA22}(C e^{i phi2} + D e^{-i phi2})]
   = Delta * f_a *)
(* Divide by v_g: *)
(* Gamma/g * [...] = Delta * etaA / g * v_g ... *)
(* Actually: f_a = (v_g/g) etaA, so Delta f_a = Delta*(v_g/g)*etaA *)
(* And the LHS: g * [...] = g * [...]. So: *)
(* (g^2/v_g)/g * [...] = (Gamma/g) * [...] = Delta * etaA *)
(* So: Gamma * [...] = Delta * etaA *)

eq13 = Gamma * (
    Exp[I thA11] (src1 + r1)                    (* WG1 x=0 *)
  + Exp[I thA12] (A Exp[I phi1] + B Exp[-I phi1])  (* WG1 x=l1 *)
  + Exp[I thA21] (src3 + r2)                    (* WG2 x=0 *)
  + Exp[I thA22] (C Exp[I phi2] + D Exp[-I phi2])  (* WG2 x=l2 *)
) - Delta * etaA;

(* Atom b equation: *)
eq14 = Gamma * (
    Exp[I thB21] (src3 + r2)                    (* WG2 x=0 *)
  + Exp[I thB22] (C Exp[I phi2] + D Exp[-I phi2])  (* WG2 x=l2 *)
  + Exp[I thB31] (0 + r3)                       (* WG3 x=0 *)
  + Exp[I thB32] (E Exp[I phi3] + F Exp[-I phi3])  (* WG3 x=l3 *)
) - Delta * etaB;

(* === SOLVE FOR WG1 INCIDENCE (src1=1, src3=0) === *)
eqsWG1 = {eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14} 
         /. {src1 -> 1, src3 -> 0};
vars = {A, B, C, D, E, F, t1, r1, t2, r2, t3, r3, etaA, etaB};

Print["Solving WG1 incidence (14 equations, 14 unknowns)..."];
solWG1 = Solve[eqsWG1 == Table[0, {14}], vars];

If[Length[solWG1] > 0,
  Print["Solution found."];
  sol = solWG1[[1]];
  Print["t1 = ", Simplify[t1 /. sol]];
  Print["r1 = ", Simplify[r1 /. sol]];
  Print["t2 = ", Simplify[t2 /. sol]];
  Print["r2 = ", Simplify[r2 /. sol]];
  Print["t3 = ", Simplify[t3 /. sol]];
  Print["r3 = ", Simplify[r3 /. sol]];
  ,
  Print["No solution found for WG1."];
];

(* === SOLVE FOR WG2 INCIDENCE (src1=0, src3=1) === *)
Print[""];
Print["Solving WG2 incidence..."];
eqsWG2 = {eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14} 
         /. {src1 -> 0, src3 -> 1};
solWG2 = Solve[eqsWG2 == Table[0, {14}], vars];

If[Length[solWG2] > 0,
  Print["Solution found."];
  sol2 = solWG2[[1]];
  (* Rename output variables with tilde *)
  Print["tt1 = ", Simplify[t1 /. sol2]];
  Print["tr1 = ", Simplify[r1 /. sol2]];
  Print["tt2 = ", Simplify[t2 /. sol2]];
  Print["tr2 = ", Simplify[r2 /. sol2]];
  Print["tt3 = ", Simplify[t3 /. sol2]];
  Print["tr3 = ", Simplify[r3 /. sol2]];
  ,
  Print["No solution found for WG2."];
];
