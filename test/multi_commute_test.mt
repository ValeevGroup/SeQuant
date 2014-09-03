(* Mathematica Test File *)

(* Test Function for commutator using function wick in MultiConfiguration *)

(* commute between two rank one excitation operator *)
Test[
    SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
    wickopts = {
       fullContract -> False,
       noCoincidences -> False,
       doSums -> True,
       doReindex -> True
       };
    p = createParticleIndex["p", occ];
    q = createParticleIndex["q", occ];
    r = createParticleIndex["r", occ];
    s = createParticleIndex["s", occ];
    t = createParticleIndex["t", occ];
    u = createParticleIndex["u", occ];
    ML1 = createSQS[{p}, {q}, inorder];
    ML2 = createSQS[{r}, {s}, inorder];
    wick[commute[ML1, ML2], {p, q, r, s}, wickopts]
    ,
    -deltaIndex[{particleIndex["s", 
    particleSpace[occupied]]}, {particleIndex["p", 
    particleSpace[occupied]]}] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] + 
    deltaIndex[{particleIndex["q", 
    particleSpace[occupied]]}, {particleIndex["r", 
    particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["s", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]] + 
    deltaIndex[{particleIndex["q", 
    particleSpace[occupied]]}, {particleIndex["r", 
    particleSpace[occupied]]}] SQS[
    particleIndex["p", particleSpace[occupied], indexType[cre]], 
    particleIndex["s", particleSpace[occupied], indexType[ann]]] - 
    deltaIndex[{particleIndex["s", 
    particleSpace[occupied]]}, {particleIndex["p", 
    particleSpace[occupied]]}] SQS[
    particleIndex["r", particleSpace[occupied], indexType[cre]], 
    particleIndex["q", particleSpace[occupied], indexType[ann]]]
    ,
    TestID->"SeQuantTest_Multi_commute_1_1"
]

(* commute between three rank one excitation operator *)
Test[
    SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
    wickopts = {
       fullContract -> False,
       noCoincidences -> False,
       doSums -> True,
       doReindex -> True
       };
    p = createParticleIndex["p", occ];
    q = createParticleIndex["q", occ];
    r = createParticleIndex["r", occ];
    s = createParticleIndex["s", occ];
    t = createParticleIndex["t", occ];
    u = createParticleIndex["u", occ];
    ML1 = createSQS[{p}, {q}, inorder];
    ML2 = createSQS[{r}, {s}, inorder];
    ML3 = createSQS[{t}, {u}, inorder];
    wick[commute[commute[ML1, ML2], ML3], {p, q, r, s, t, u}, wickopts]
    ,
    deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["u", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["t", particleSpace[occupied], indexType[ket]]] - 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["u", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["s", particleSpace[occupied], indexType[bra]], 
    particleIndex["t", particleSpace[occupied], indexType[ket]]] + 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["t", 
     particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]] - 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["t", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] + 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["t", 
     particleSpace[occupied]]}] SQS[
    particleIndex["p", particleSpace[occupied], indexType[cre]], 
    particleIndex["u", particleSpace[occupied], indexType[ann]]] - 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["t", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] SQS[
    particleIndex["r", particleSpace[occupied], indexType[cre]], 
    particleIndex["u", particleSpace[occupied], indexType[ann]]] + 
    deltaIndex[{particleIndex["s", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["u", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] SQS[
    particleIndex["t", particleSpace[occupied], indexType[cre]], 
    particleIndex["q", particleSpace[occupied], indexType[ann]]] - 
    deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] deltaIndex[{particleIndex["u", 
     particleSpace[occupied]]}, {particleIndex["p", 
     particleSpace[occupied]]}] SQS[
    particleIndex["t", particleSpace[occupied], indexType[cre]], 
    particleIndex["s", particleSpace[occupied], indexType[ann]]]
    ,
    TestID->"SeQuantTest_Multi_commute_1_1_1"
]