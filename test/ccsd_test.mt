(* Mathematica Test File *)


(* Test CCSD Energy equations *)
Test[
    i = createParticleIndex["i", occ];
    j = createParticleIndex["j", occ];
    k = createParticleIndex["k", occ];
    l = createParticleIndex["l", occ];
    m = createParticleIndex["m", occ];
    n = createParticleIndex["n", occ];
    Subscript[i, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", occ];
    Subscript[j, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", occ];
    Subscript[k, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(1\)]\)", occ];
    Subscript[l, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(1\)]\)", occ];
    Subscript[i, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", occ];
    Subscript[j, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(2\)]\)", occ];
    Subscript[k, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(2\)]\)", occ];
    Subscript[l, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(2\)]\)", occ];
    Subscript[i, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(3\)]\)", occ];
    Subscript[j, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(3\)]\)", occ];
    Subscript[k, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(3\)]\)", occ];
    Subscript[l, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(3\)]\)", occ];
    Subscript[i, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(4\)]\)", occ];
    Subscript[j, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(4\)]\)", occ];
    Subscript[k, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(4\)]\)", occ];
    Subscript[l, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(4\)]\)", occ];
    a = createParticleIndex["a", virt];
    b = createParticleIndex["b", virt];
    c = createParticleIndex["c", virt];
    d = createParticleIndex["d", virt];
    e = createParticleIndex["e", virt];
    f = createParticleIndex["f", virt];
    Subscript[a, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", virt];
    Subscript[b, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", virt];
    Subscript[c, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(1\)]\)", virt];
    Subscript[d, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(1\)]\)", virt];
    Subscript[a, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", virt];
    Subscript[b, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(2\)]\)", virt];
    Subscript[c, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(2\)]\)", virt];
    Subscript[d, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(2\)]\)", virt];
    Subscript[a, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(3\)]\)", virt];
    Subscript[b, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(3\)]\)", virt];
    Subscript[c, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(3\)]\)", virt];
    Subscript[d, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(3\)]\)", virt];
    Subscript[a, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(4\)]\)", virt];
    Subscript[b, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(4\)]\)", virt];
    Subscript[c, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(4\)]\)", virt];
    Subscript[d, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(4\)]\)", virt];
    \[Kappa] = createParticleIndex["\[Kappa]", allany];
    \[Lambda] = createParticleIndex["\[Lambda]", allany];
    \[Mu] = createParticleIndex["\[Mu]", allany];
    \[Nu] = createParticleIndex["\[Nu]", allany];
    g2 = createSQM["\!\(\*OverscriptBox[\(g\), \(_\)]\)", {\[Mu], \[Nu]}, {\[Kappa], \[Lambda]}, antisymm];
    ga2 = createSQS[{\[Mu], \[Nu]}, {\[Kappa], \[Lambda]}];
    H2 = 1/4 g2*ga2;
    f1 = createSQM["F", {\[Mu]}, {\[Nu]}, antisymm];
    ga1 = createSQS[{\[Mu]}, {\[Nu]}];
    H1 = f1*ga1;
    H = H1 + H2;
    wikopt = {fullContract -> False, noCoincidences -> False, 
       useDensity -> True, doSums -> True, doReindex -> True};
    Subscript[T1, 1] = 
      createSQM["t", {Subscript[a, 1]}, {Subscript[i, 1]}, antisymm]*
       createSQS[{Subscript[a, 1]}, {Subscript[i, 1]}];
    Subscript[T1, 2] = 
      createSQM["t", {Subscript[a, 2]}, {Subscript[i, 2]}, antisymm]*
       createSQS[{Subscript[a, 2]}, {Subscript[i, 2]}];
    Subscript[T1, 3] = 
      createSQM["t", {Subscript[a, 3]}, {Subscript[i, 3]}, antisymm]*
       createSQS[{Subscript[a, 3]}, {Subscript[i, 3]}];
    Subscript[T1, 4] = 
      createSQM["t", {Subscript[a, 4]}, {Subscript[i, 4]}, antisymm]*
       createSQS[{Subscript[a, 4]}, {Subscript[i, 4]}];
    Subscript[T2, 1] = 
      1/4*createSQM[
        "t", {Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}, antisymm]*
       createSQS[{Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}];
    Subscript[T2, 2] = 
      1/4*createSQM[
        "t", {Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}, antisymm]*
       createSQS[{Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}];
    Subscript[T2, 3] = 
      1/4*createSQM[
        "t", {Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}, antisymm]*
       createSQS[{Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}];
    Subscript[T2, 4] = 
      1/4*createSQM[
        "t", {Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}, antisymm]*
       createSQS[{Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}];

    (* projection operators *)
    P1 = createSQS[{i}, {a}];
    P2 = createSQS[{i, j}, {a, b}];

    (* The function processes expressions after Wick's theorem *)
    process[expr_, extIndices_List: {}, 
       params_List: defaultTheoryParams] :=
        Module[ {result},
            result = expr;
            If[ useBrillouin /. params,
                result = brillouin[result]
            ];
            If[ !keepDisconnected /. params,
                result = removeDisconnectedTerms[result, extIndices]
            ];
            If[ doSimplify /. params,
                result = Simplify[result]
            ];
            Return[result];
        ];
    Subscript[T, 1] = Subscript[T1, 1] + Subscript[T2, 1];
    Subscript[T, 2] = Subscript[T1, 2] + Subscript[T2, 2];
    Subscript[T, 3] = Subscript[T1, 3] + Subscript[T2, 3];
    Subscript[T, 4] = Subscript[T1, 4] + Subscript[T2, 4];
    defaultTheoryParams =
      {
       useBrillouin -> False,
       keepDisconnected -> False,
       doSimplify -> True
       };
    theoryParams = defaultTheoryParams;
    SeQuantDebugLevel = 0;
    tmp1 = wick[H ** Subscript[T, 1], {}];
    tmp1 = process[tmp1, {}, theoryParams];
    
    (* 1/2 H**T**T->0 *)
    tmp2 = wick[H ** Subscript[T, 1] ** Subscript[T, 2], {}];
    tmp2 /= 2;
    tmp2 = process[tmp2, {}, theoryParams];
    tmp = tmp1 + tmp2
    ,
    SQM[OHead["F", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    1/4 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] +
    1/2 SQM[OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
    particleSpace[virtual], indexType[ket]]]
    ,
    TestID->"SeQuantTest_ccsd_energy"
]


(* Test CCSD T1 Amplitude Equations *)
Test[
    ClearAll["Global`*"];
    i = createParticleIndex["i", occ];
    j = createParticleIndex["j", occ];
    k = createParticleIndex["k", occ];
    l = createParticleIndex["l", occ];
    m = createParticleIndex["m", occ];
    n = createParticleIndex["n", occ];
    Subscript[i, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", occ];
    Subscript[j, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", occ];
    Subscript[k, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(1\)]\)", occ];
    Subscript[l, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(1\)]\)", occ];
    Subscript[i, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", occ];
    Subscript[j, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(2\)]\)", occ];
    Subscript[k, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(2\)]\)", occ];
    Subscript[l, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(2\)]\)", occ];
    Subscript[i, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(3\)]\)", occ];
    Subscript[j, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(3\)]\)", occ];
    Subscript[k, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(3\)]\)", occ];
    Subscript[l, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(3\)]\)", occ];
    Subscript[i, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(4\)]\)", occ];
    Subscript[j, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(4\)]\)", occ];
    Subscript[k, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(4\)]\)", occ];
    Subscript[l, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(4\)]\)", occ];
    a = createParticleIndex["a", virt];
    b = createParticleIndex["b", virt];
    c = createParticleIndex["c", virt];
    d = createParticleIndex["d", virt];
    e = createParticleIndex["e", virt];
    f = createParticleIndex["f", virt];
    Subscript[a, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", virt];
    Subscript[b, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", virt];
    Subscript[c, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(1\)]\)", virt];
    Subscript[d, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(1\)]\)", virt];
    Subscript[a, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", virt];
    Subscript[b, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(2\)]\)", virt];
    Subscript[c, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(2\)]\)", virt];
    Subscript[d, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(2\)]\)", virt];
    Subscript[a, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(3\)]\)", virt];
    Subscript[b, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(3\)]\)", virt];
    Subscript[c, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(3\)]\)", virt];
    Subscript[d, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(3\)]\)", virt];
    Subscript[a, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(4\)]\)", virt];
    Subscript[b, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(4\)]\)", virt];
    Subscript[c, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(4\)]\)", virt];
    Subscript[d, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(4\)]\)", virt];
    \[Kappa] = createParticleIndex["\[Kappa]", allany];
    \[Lambda] = createParticleIndex["\[Lambda]", allany];
    \[Mu] = createParticleIndex["\[Mu]", allany];
    \[Nu] = createParticleIndex["\[Nu]", allany];

    (* Hamiltonian *)
    g2 = createSQM[
       "\!\(\*OverscriptBox[\(g\), \(_\)]\)", {\[Mu], \[Nu]}, {\[Kappa], \
\[Lambda]}, antisymm];
    ga2 = createSQS[{\[Mu], \[Nu]}, {\[Kappa], \[Lambda]}];
    H2 = 1/4 g2*ga2;
    f1 = createSQM["F", {\[Mu]}, {\[Nu]}, antisymm];
    ga1 = createSQS[{\[Mu]}, {\[Nu]}];
    H1 = f1*ga1;
    H = H1 + H2;
    wikopt = {fullContract -> False, noCoincidences -> False, 
       useDensity -> True, doSums -> True, doReindex -> True};

    (*
    Cluster operators
    *)
    Subscript[T1, 1] = 
      createSQM["t", {Subscript[a, 1]}, {Subscript[i, 1]}, antisymm]*
       createSQS[{Subscript[a, 1]}, {Subscript[i, 1]}];
    Subscript[T1, 2] = 
      createSQM["t", {Subscript[a, 2]}, {Subscript[i, 2]}, antisymm]*
       createSQS[{Subscript[a, 2]}, {Subscript[i, 2]}];
    Subscript[T1, 3] = 
      createSQM["t", {Subscript[a, 3]}, {Subscript[i, 3]}, antisymm]*
       createSQS[{Subscript[a, 3]}, {Subscript[i, 3]}];
    Subscript[T1, 4] = 
      createSQM["t", {Subscript[a, 4]}, {Subscript[i, 4]}, antisymm]*
       createSQS[{Subscript[a, 4]}, {Subscript[i, 4]}];
    Subscript[T2, 1] = 
      1/4*createSQM[
        "t", {Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}, antisymm]*
       createSQS[{Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}];
    Subscript[T2, 2] = 
      1/4*createSQM[
        "t", {Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}, antisymm]*
       createSQS[{Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}];
    Subscript[T2, 3] = 
      1/4*createSQM[
        "t", {Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}, antisymm]*
       createSQS[{Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}];
    Subscript[T2, 4] = 
      1/4*createSQM[
        "t", {Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}, antisymm]*
       createSQS[{Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}];

    (* projection operators *)
    P1 = createSQS[{i}, {a}];
    P2 = createSQS[{i, j}, {a, b}];

    (* The function processes expressions after Wick's theorem *)
    process[expr_, extIndices_List: {}, 
       params_List: defaultTheoryParams] :=
        Module[ {result},
            result = expr;
            If[ useBrillouin /. params,
                result = brillouin[result]
            ];
            If[ !keepDisconnected /. params,
                result = removeDisconnectedTerms[result, extIndices]
            ];
            If[ doSimplify /. params,
                result = Simplify[result]
            ];
            Return[result];
        ];
    Subscript[T, 1] = Subscript[T1, 1] + Subscript[T2, 1];
    Subscript[T, 2] = Subscript[T1, 2] + Subscript[T2, 2];
    Subscript[T, 3] = Subscript[T1, 3] + Subscript[T2, 3];
    Subscript[T, 4] = Subscript[T1, 4] + Subscript[T2, 4];
    defaultTheoryParams =
      {
       useBrillouin -> False,
       keepDisconnected -> False,
       doSimplify -> True
       };
    theoryParams = defaultTheoryParams;
    SeQuantDebugLevel = 0;
    T1ExtIndices = {i, a};
    tmp1 = wick[P1 ** H ** Subscript[T, 1], T1ExtIndices];
    tmp1 = process[tmp1, T1ExtIndices, theoryParams];
    tmp2 = wick[P1 ** H ** Subscript[T, 1] ** Subscript[T, 2], T1ExtIndices];
    tmp2 = process[tmp2, T1ExtIndices, theoryParams];
    tmp2 /= 2;
    tmp3 = wick[P1 ** H ** Subscript[T1, 1] ** Subscript[T1, 2] ** Subscript[T1, 3], T1ExtIndices];
    tmp3 = process[tmp3, T1ExtIndices, theoryParams];
    tmp3 /= 6;
    tmp1 + tmp2 + tmp3
    ,
    -SQM[OHead["F", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["F", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]] SQM[
     OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]]] - 
    SQM[OHead["t", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
      particleSpace[occupied], indexType[ket]]] SQM[
     OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
      particleSpace[virtual], indexType[ket]]] + 
    SQM[OHead["F", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]] - 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]] + 
    1/2 (-2 SQM[OHead["F", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[ket]]] SQM[
        OHead["t", indexSymm[-1]], 
        particleIndex["a", particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[ket]]] SQM[
        OHead["t", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["i", particleSpace[occupied], indexType[ket]]] - 
      2 SQM[OHead["t", indexSymm[-1]], 
        particleIndex["a", particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[ket]]] SQM[
        OHead["t", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[ket]]] SQM[
        OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["i", particleSpace[occupied], indexType[ket]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[ket]]] + 
      2 SQM[OHead["t", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[ket]]] SQM[
        OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[ket]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
         particleSpace[virtual], indexType[ket]]] SQM[
        OHead["t", indexSymm[-1]], 
        particleIndex["a", particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["i", particleSpace[occupied], indexType[ket]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[ket]]] + 
      SQM[OHead["t", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["i", particleSpace[occupied], 
         indexType[
          ket]]] (2 SQM[OHead["t", indexSymm[-1]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
            particleSpace[virtual], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
            particleSpace[occupied], indexType[ket]]] SQM[
           OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
           particleIndex["a", particleSpace[virtual], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
            particleSpace[occupied], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
            particleSpace[virtual], indexType[ket]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
            particleSpace[virtual], indexType[ket]]] - 
         SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
            indexSymm[-1]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
            particleSpace[occupied], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
            particleSpace[occupied], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
            particleSpace[virtual], indexType[ket]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
            particleSpace[virtual], indexType[ket]]] SQM[
           OHead["t", indexSymm[-1]], 
           particleIndex["a", particleSpace[virtual], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
            particleSpace[virtual], indexType[bra]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
            particleSpace[occupied], indexType[ket]], 
           particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
            particleSpace[occupied], indexType[ket]]]) - 
      SQM[OHead["t", indexSymm[-1]], 
        particleIndex["a", particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[ket]]] SQM[
        OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[ket]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
         particleSpace[virtual], indexType[ket]]] SQM[
        OHead["t", indexSymm[-1]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
         particleSpace[virtual], indexType[bra]], 
        particleIndex["i", particleSpace[occupied], indexType[ket]], 
        particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
         particleSpace[occupied], indexType[ket]]]) + 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
     particleIndex["a", particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
      particleSpace[virtual], indexType[ket]]] SQM[
     OHead["t", indexSymm[-1]], 
     particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
      particleSpace[virtual], indexType[bra]], 
     particleIndex["i", particleSpace[occupied], indexType[ket]], 
     particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
      particleSpace[occupied], indexType[ket]]]
    ,
    TestID->"SeQuantTest_ccsd_T1amplitude"
];

(* Test CCSD T2 Amplitude Equations *)
Test[
    ClearAll["Global`*"];
    i = createParticleIndex["i", occ];
    j = createParticleIndex["j", occ];
    k = createParticleIndex["k", occ];
    l = createParticleIndex["l", occ];
    m = createParticleIndex["m", occ];
    n = createParticleIndex["n", occ];
    Subscript[i, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", occ];
    Subscript[j, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", occ];
    Subscript[k, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(1\)]\)", occ];
    Subscript[l, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(1\)]\)", occ];
    Subscript[i, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", occ];
    Subscript[j, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(2\)]\)", occ];
    Subscript[k, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(2\)]\)", occ];
    Subscript[l, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(2\)]\)", occ];
    Subscript[i, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(3\)]\)", occ];
    Subscript[j, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(3\)]\)", occ];
    Subscript[k, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(3\)]\)", occ];
    Subscript[l, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(3\)]\)", occ];
    Subscript[i, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(i\), \(4\)]\)", occ];
    Subscript[j, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(j\), \(4\)]\)", occ];
    Subscript[k, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(k\), \(4\)]\)", occ];
    Subscript[l, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(l\), \(4\)]\)", occ];
    a = createParticleIndex["a", virt];
    b = createParticleIndex["b", virt];
    c = createParticleIndex["c", virt];
    d = createParticleIndex["d", virt];
    e = createParticleIndex["e", virt];
    f = createParticleIndex["f", virt];
    Subscript[a, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", virt];
    Subscript[b, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", virt];
    Subscript[c, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(1\)]\)", virt];
    Subscript[d, 1] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(1\)]\)", virt];
    Subscript[a, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", virt];
    Subscript[b, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(2\)]\)", virt];
    Subscript[c, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(2\)]\)", virt];
    Subscript[d, 2] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(2\)]\)", virt];
    Subscript[a, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(3\)]\)", virt];
    Subscript[b, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(3\)]\)", virt];
    Subscript[c, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(3\)]\)", virt];
    Subscript[d, 3] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(3\)]\)", virt];
    Subscript[a, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(a\), \(4\)]\)", virt];
    Subscript[b, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(b\), \(4\)]\)", virt];
    Subscript[c, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(c\), \(4\)]\)", virt];
    Subscript[d, 4] = 
      createParticleIndex["\!\(\*SubscriptBox[\(d\), \(4\)]\)", virt];
    \[Kappa] = createParticleIndex["\[Kappa]", allany];
    \[Lambda] = createParticleIndex["\[Lambda]", allany];
    \[Mu] = createParticleIndex["\[Mu]", allany];
    \[Nu] = createParticleIndex["\[Nu]", allany];

    (* Hamiltonian *)
    g2 = createSQM[
       "\!\(\*OverscriptBox[\(g\), \(_\)]\)", {\[Mu], \[Nu]}, {\[Kappa], \
\[Lambda]}, antisymm];
    ga2 = createSQS[{\[Mu], \[Nu]}, {\[Kappa], \[Lambda]}];
    H2 = 1/4 g2*ga2;
    f1 = createSQM["F", {\[Mu]}, {\[Nu]}, antisymm];
    ga1 = createSQS[{\[Mu]}, {\[Nu]}];
    H1 = f1*ga1;
    H = H1 + H2;
    wikopt = {fullContract -> False, noCoincidences -> False, 
       useDensity -> True, doSums -> True, doReindex -> True};

    (*
    Cluster operators
    *)
    Subscript[T1, 1] = 
      createSQM["t", {Subscript[a, 1]}, {Subscript[i, 1]}, antisymm]*
       createSQS[{Subscript[a, 1]}, {Subscript[i, 1]}];
    Subscript[T1, 2] = 
      createSQM["t", {Subscript[a, 2]}, {Subscript[i, 2]}, antisymm]*
       createSQS[{Subscript[a, 2]}, {Subscript[i, 2]}];
    Subscript[T1, 3] = 
      createSQM["t", {Subscript[a, 3]}, {Subscript[i, 3]}, antisymm]*
       createSQS[{Subscript[a, 3]}, {Subscript[i, 3]}];
    Subscript[T1, 4] = 
      createSQM["t", {Subscript[a, 4]}, {Subscript[i, 4]}, antisymm]*
       createSQS[{Subscript[a, 4]}, {Subscript[i, 4]}];
    Subscript[T2, 1] = 
      1/4*createSQM[
        "t", {Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}, antisymm]*
       createSQS[{Subscript[a, 1], Subscript[b, 1]}, {Subscript[i, 1], 
         Subscript[j, 1]}];
    Subscript[T2, 2] = 
      1/4*createSQM[
        "t", {Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}, antisymm]*
       createSQS[{Subscript[a, 2], Subscript[b, 2]}, {Subscript[i, 2], 
         Subscript[j, 2]}];
    Subscript[T2, 3] = 
      1/4*createSQM[
        "t", {Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}, antisymm]*
       createSQS[{Subscript[a, 3], Subscript[b, 3]}, {Subscript[i, 3], 
         Subscript[j, 3]}];
    Subscript[T2, 4] = 
      1/4*createSQM[
        "t", {Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}, antisymm]*
       createSQS[{Subscript[a, 4], Subscript[b, 4]}, {Subscript[i, 4], 
         Subscript[j, 4]}];

    (* projection operators *)
    P1 = createSQS[{i}, {a}];
    P2 = createSQS[{i, j}, {a, b}];

    (* The function processes expressions after Wick's theorem *)
    process[expr_, extIndices_List: {}, 
       params_List: defaultTheoryParams] :=
        Module[ {result},
            result = expr;
            If[ useBrillouin /. params,
                result = brillouin[result]
            ];
            If[ !keepDisconnected /. params,
                result = removeDisconnectedTerms[result, extIndices]
            ];
            If[ doSimplify /. params,
                result = Simplify[result]
            ];
            Return[result];
        ];
    Subscript[T, 1] = Subscript[T1, 1] + Subscript[T2, 1];
    Subscript[T, 2] = Subscript[T1, 2] + Subscript[T2, 2];
    Subscript[T, 3] = Subscript[T1, 3] + Subscript[T2, 3];
    Subscript[T, 4] = Subscript[T1, 4] + Subscript[T2, 4];
    defaultTheoryParams =
      {
       useBrillouin -> False,
       keepDisconnected -> False,
       doSimplify -> True
       };
    theoryParams = defaultTheoryParams;
    SeQuantDebugLevel = 0;
    T2ExtIndices = {i, j, a, b};
    tmp1 = wick[P2 ** H ** Subscript[T, 1], T2ExtIndices];
    tmp1 = process[tmp1, T2ExtIndices, theoryParams];
    tmp2 = wick[P2 ** H ** Subscript[T, 1] ** Subscript[T, 2], T2ExtIndices];
    tmp2 = process[tmp2, T2ExtIndices, theoryParams];
    tmp2 /= 2;
    tmp3 = wick[P2 ** H ** Subscript[T1, 1] ** Subscript[T1, 2] ** Subscript[T1, 3], T2ExtIndices];
    tmp3 = process[tmp3, T2ExtIndices, theoryParams];
    tmp3 /= 6;
    tmp4 = wick[P2 ** H ** Subscript[T2, 1] ** Subscript[T1, 2] ** Subscript[T1, 3], T2ExtIndices];
    tmp4 = process[tmp4, T2ExtIndices, theoryParams];
    tmp4 /= 2;
    tmp5 = wick[P2 ** H ** Subscript[T1, 1] ** Subscript[T1, 2] ** Subscript[T1, 3] ** Subscript[T1, 4], T2ExtIndices];
    tmp5 = process[tmp5, T2ExtIndices, theoryParams];
    tmp5 /= 24;
    tmp1 + tmp2 + tmp3 + tmp4 + tmp5
    ,
    SQM[OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] - 
    SQM[OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] - 
    SQM[OHead["t", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    1/6 (-6 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] + 
    6 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], 
       indexType[
        ket]]] (SQM[OHead["t", indexSymm[-1]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], 
          indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] + 
       SQM[OHead["t", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], 
          indexType[
           ket]]] (SQM[OHead["t", indexSymm[-1]], 
            particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
             particleSpace[virtual], indexType[bra]], 
            particleIndex["j", particleSpace[occupied], 
             indexType[ket]]] SQM[
            OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
             indexSymm[-1]], 
            particleIndex["b", particleSpace[virtual], 
             indexType[bra]], 
            
            particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
             particleSpace[occupied], indexType[bra]], 
            particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
             particleSpace[virtual], indexType[ket]], 
            particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
             particleSpace[virtual], indexType[ket]]] - 
          SQM[OHead["t", indexSymm[-1]], 
            particleIndex["b", particleSpace[virtual], 
             indexType[bra]], 
            particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
             particleSpace[occupied], indexType[ket]]] SQM[
            OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
             indexSymm[-1]], 
            particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
             particleSpace[occupied], indexType[bra]], 
            particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
             particleSpace[occupied], indexType[bra]], 
            particleIndex["j", particleSpace[occupied], 
             indexType[ket]], 
            particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
             particleSpace[virtual], indexType[ket]]]))) + 
    SQM[OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
    particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
    particleSpace[virtual], indexType[ket]]] - 
    SQM[OHead["F", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["F", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(j\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["F", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["F", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[ket]]] + 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
    particleSpace[occupied], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
    particleSpace[virtual], 
    indexType[
     ket]]] (-2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[
        ket]]] (2 SQM[OHead["t", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
          particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] + 
       SQM[OHead["t", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
          particleSpace[virtual], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], 
          indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[ket]]]) - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]]) + 
    1/2 (-SQM[OHead["t", indexSymm[-1]], 
       particleIndex["a", particleSpace[virtual], indexType[bra]], 
       particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
        particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[
        ket]]] (-SQM[OHead["t", indexSymm[-1]], 
          particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
           particleSpace[virtual], indexType[bra]], 
          particleIndex["i", particleSpace[occupied], 
           indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
          particleSpace[virtual], indexType[ket]]] - 
       2 SQM[OHead["t", indexSymm[-1]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] + 
       2 SQM[OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] - 
       2 SQM[OHead["F", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] + 
       SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
          indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["i", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[ket]]]) + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[
        ket]]] (SQM[OHead["t", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
          particleSpace[virtual], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], 
          indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
          particleSpace[virtual], indexType[ket]]] + 
       2 SQM[OHead["t", indexSymm[-1]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] - 
       2 SQM[OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] SQM[
         OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] + 
       2 SQM[OHead["F", indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]]] - 
       SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
          indexSymm[-1]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[bra]], 
         particleIndex["j", particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
          particleSpace[virtual], indexType[ket]]] SQM[
         OHead["t", indexSymm[-1]], 
         particleIndex["a", particleSpace[virtual], indexType[bra]], 
         particleIndex["b", particleSpace[virtual], indexType[bra]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
          particleSpace[occupied], indexType[ket]], 
         particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
          particleSpace[occupied], indexType[ket]]]) - 
    2 SQM[OHead["F", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["F", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], 
       indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]],
       particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] + 
    2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]],
       particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], 
       indexType[ket]]] SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", 
       indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]]] + 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]] - 
    SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["a", particleSpace[virtual], indexType[bra]], 
      particleIndex["b", particleSpace[virtual], indexType[bra]], 
      particleIndex["i", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(1\)]\)", 
       particleSpace[occupied], indexType[ket]]] SQM[
      OHead["t", indexSymm[-1]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["\!\(\*SubscriptBox[\(a\), \(2\)]\)", 
       particleSpace[virtual], indexType[bra]], 
      particleIndex["j", particleSpace[occupied], indexType[ket]], 
      particleIndex["\!\(\*SubscriptBox[\(i\), \(2\)]\)", 
       particleSpace[occupied], indexType[ket]]]) + 
    1/2 SQM[OHead["\!\(\*OverscriptBox[\(g\), \(_\)]\)", indexSymm[-1]], 
    particleIndex["a", particleSpace[virtual], indexType[bra]], 
    particleIndex["b", particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]], 
    particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
    particleSpace[virtual], indexType[ket]]] SQM[
    OHead["t", indexSymm[-1]], 
    particleIndex["\!\(\*SubscriptBox[\(a\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["\!\(\*SubscriptBox[\(b\), \(1\)]\)", 
    particleSpace[virtual], indexType[bra]], 
    particleIndex["i", particleSpace[occupied], indexType[ket]], 
    particleIndex["j", particleSpace[occupied], indexType[ket]]]
    ,
    TestID->"SeQuantTest_ccsd_T2amplitude"
];