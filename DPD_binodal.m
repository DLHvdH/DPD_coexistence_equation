(* 
    Module:	DPD_binodal.m
    Authors:	Dingeman L.H. van der Haven, Stephan Koehler, Eduard Schreiner,
		and Pieter J. in 't Veld
    Date:	March 3, 2021
    Purpose:	Binodal curve equations
*)

Module[{
  at, nt, B2,
  kappa, lambda, delta,
  aref = 25,
  c, cerr, k = {
    {{0.948727, 0.0079}, {0.4127699, 0.0052}},
    {{-0.2666, 0.0410}, {1.443, 0.111}, {-0.6669, 0.0871}, {0.08599, 
      0.02057}, {1.869, 0.036}, {-2.339, 0.143}},
    {{0.04474, 0.00262}, {1.740, 0.043}, {0.8931, 0.0346}, {1.287, 
      0.049}, {5.924, 1.922}, {1.127, 0.054}, {-1.519, 
      0.273}, {0.4926, 0.0384}},
    {{1.141, 0.063}, {-1.077, 0.170}, {0.7211, 0.1338}, {-0.09151, 
      0.03161}, {-1.006, 0.027}, {1.758, 0.113}},
    {{2.33124, 0.164}, {0.130069, 0.0101}, {0.00263761, 0.000145}, { 
      1.44442, 0.013 }, {0.362353, 0.0422}, {-0.582265, 0.0620}}}
  },
 c = Map[Map[#[[1]] &, #] &, k];
 cerr = Map[Map[#[[2]] &, #] &, k];
 
 at[aii_, ajj_] := (aii + ajj)/(2 aref);
 
 nt[ni_, nj_] := 1 - 2 Sqrt[ni nj]/(ni + nj);
 
 B2[a_] = Integrate[2*Pi*(1 - Exp[-a/2*(1 - r)^2])*r^2, {r, 0, 1}];
 
 phicVDH[aii_, ajj_, ni_, nj_] := 
  1/(1 + (B2[aii]/B2[ajj])^c[[1, 1]]*(ni/nj)^c[[1, 2]]);

 rhoiiVDH[aii_] := 12.4965358024165/aii^0.438805064915560; (* at p = 24.7 *)
 
 kappa[aii_, ajj_, ni_, nj_] :=
  
  Sum[c[[2, i]] at[aii, ajj]^(i - 1), {i, 4}]*(1 + 
     Sum[c[[2, i + 4]] (2 phicVDH[aii, ajj, 1, 1] - 1)^i, {i, 2}])*
   (1 + c[[3, 1]] Log[nj]^c[[3, 2]])/(1 + c[[3, 3]] Log[ni]^c[[3, 4]])*
   (1 + c[[3, 5]] (Log[ni] Log[nj])^c[[3, 6]])*
   Exp[c[[3, 7]] (Log[ni] Log[nj])^c[[3, 8]]];
 
 lambda[aii_, ajj_, ni_, nj_] :=
  
  1 + Sum[c[[4, i]] at[aii, ajj]^(i - 1), {i, 
      4}]*(1 + 
       Sum[c[[4, i + 4]] (2 phicVDH[aii, ajj, 1, 1] - 1)^i, {i, 
         2}])/(ni^(nj/3) + nj^(ni/9) - 1);
 
 delta[aii_, ajj_, ni_, nj_] := Block[{dm, am = Sqrt[aii ajj]},
   dm = am + Sum[c[[5, i]] am^(i - 1), {i, 3}];
   am*(1 + (dm - 
          am)/(2 am)*(ni^(c[[5, 4]]/2) + nj^(c[[5, 4]]/2))^2/(2 ni^
           c[[5, 4]] nj^c[[5, 4]])*
       (1 + Sum[c[[5, i + 4]] nt[ni, nj]^i, {i, 2}]))];
 
 aijVDH[phi_, aii_, ajj_, ni_, nj_] := 
  kappa[aii, ajj, ni, nj] Log[1/phi]^lambda[aii, ajj, ni, nj] +
   kappa[ajj, aii, nj, ni] Log[1/(1 - phi)]^
     lambda[ajj, aii, nj, ni] +
   delta[aii, ajj, ni, nj]
 ]
