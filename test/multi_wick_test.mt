(* Mathematica Test File *)

(* Test Function for wick function in MultiConfiguration *)

(*

normal contraction, using all occupied orbitals

*)


(* wick between rank one and rank one operator *)
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
	ML1 = mSQS[normalOrder[True], 
 SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  particleIndex["q", particleSpace[occupied], indexType[ann]]]];
  
  	ML2 = mSQS[normalOrder[True], 
 SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
  particleIndex["s", particleSpace[occupied], indexType[ann]]]];
  
	wick[ML1 ** ML2, {p, q, r, s}, wickopts]
	,	
	deltaIndex[particleIndex["q", 
     particleSpace[occupied]], particleIndex["r", 
     particleSpace[occupied]]] SQM[OHead["\[Lambda]", indexSymm[-1]],
    particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] - 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
   OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] + 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
  particleIndex["q", particleSpace[occupied], indexType[bra]], 
  particleIndex["s", particleSpace[occupied], indexType[bra]], 
  particleIndex["p", particleSpace[occupied], indexType[ket]], 
  particleIndex["r", particleSpace[occupied], indexType[ket]]] + 
 deltaIndex[particleIndex["q", 
     particleSpace[occupied]], particleIndex["r", 
     particleSpace[occupied]]] SQS[
   particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[ann]]] - 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]] SQS[
   particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[ann]]] - 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] SQS[
   particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]] + 
 SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  particleIndex["r", particleSpace[occupied], indexType[cre]], 
  particleIndex["s", particleSpace[occupied], indexType[ann]], 
  particleIndex["q", particleSpace[occupied], indexType[ann]]]
	,	
	TestID->"SeQuantTest_Multi_wick_1_1"
]


(* wick between rank two and rank two operator *)
Test[
    SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
    SeQuantDebugLevel = 0;
    wickopts = {
       fullContract -> False,
       noCoincidences -> False,
       doSums -> True,
       doReindex -> True
       };

    (* One-particle indices *)
    p = createParticleIndex["p", occ];
    q = createParticleIndex["q", occ];
    r = createParticleIndex["r", occ];
    s = createParticleIndex["s", occ];
    t = createParticleIndex["t", occ];
    u = createParticleIndex["u", occ];
    v = createParticleIndex["v", occ];
    w = createParticleIndex["w", occ];
    ML4 = createSQS[{p, q}, {r, s}, noorder];
    ML5 = createSQS[{t, u}, {v, w}, noorder];
    ML7 = createSQS[{p, q}, {v, w}, noorder];
    ML8 = createSQS[{p, q, t}, {s, v, w}, noorder];
    ML9 = createSQS[{p, q, t}, {r, v, w}, noorder];
    ML10 = createSQS[{p, q, t, u}, {r, s, v, w}, noorder];
    ML11 = createSQS[{p, q, u}, {s, v, w}, noorder];
    ML12 = createSQS[{p, q, u}, {r, v, w}, noorder];
    res = wick[ML4 ** ML5, {p, q, r, s, t, u, v, w}, wickopts];
    nor = - normalOrderedForm[ML7]*deltaIndex[r, u]*deltaIndex[s, t] + 
       normalOrderedForm[ML7]*deltaIndex[r, t]*deltaIndex[s, u] + 
       normalOrderedForm[ML8]*deltaIndex[r, u] - 
       normalOrderedForm[ML9]*deltaIndex[s, u] + 
       normalOrderedForm[ML10] - 
       normalOrderedForm[ML11]*deltaIndex[r, t] + 
       normalOrderedForm[ML12]*deltaIndex[s, t];
    nor = Map[Distribute, nor];
    res === nor
    ,
    True
    ,
    TestID->"SeQuantTest_Multi_wick_2_2"
]


(* wick between rank three and rank one operator *)
Test[
    SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
    SeQuantDebugLevel = 0;
    wickopts = {
       fullContract -> False,
       noCoincidences -> False,
       doSums -> True,
       doReindex -> True
       };

    (* One-particle indices *)
    p = createParticleIndex["p", occ];
    q = createParticleIndex["q", occ];
    r = createParticleIndex["r", occ];
    s = createParticleIndex["s", occ];
    t = createParticleIndex["t", occ];
    u = createParticleIndex["u", occ];
    v = createParticleIndex["v", occ];
    w = createParticleIndex["w", occ];
    ML1 = createSQS[{p, q, r}, {s, t, u}, noorder];
    ML2 = createSQS[{v}, {w}, noorder];
    ML3 = createSQS[{p, q, r}, {t, u, w}, noorder];
    ML4 = createSQS[{p, q, r}, {s, t, w}, noorder];
    ML5 = createSQS[{p, q, r}, {s, u, w}, noorder];
    ML6 = createSQS[{p, q, r, v}, {s, t, u, w}, noorder];
    res = wick[ML1 ** ML2, {p, q, r, s, t, u, v, w}, wickopts];
    nor = normalOrderedForm[ML3]*deltaIndex[s,v] + normalOrderedForm[ML4]*deltaIndex[u,v] - normalOrderedForm[ML5]*deltaIndex[t,v] + normalOrderedForm[ML6];
    nor = Map[Distribute, nor];
    res === nor
    ,
    True
    ,
    TestID->"SeQuantTest_Multi_wick_3_1"
]

(* wick between rank three and rank two operator *)
(*Test[
    SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
    SeQuantDebugLevel = 0;
    wickopts = {
       fullContract -> False,
       noCoincidences -> False,
       doSums -> True,
       doReindex -> True
       };

    (* One-particle indices *)
    p = createParticleIndex["p", occ];
    q = createParticleIndex["q", occ];
    r = createParticleIndex["r", occ];
    s = createParticleIndex["s", occ];
    t = createParticleIndex["t", occ];
    u = createParticleIndex["u", occ];
    v = createParticleIndex["v", occ];
    w = createParticleIndex["w", occ];
    a = createParticleIndex["a", occ];
	b = createParticleIndex["b", occ];
	
    ML1 = createSQS[{p, q, r}, {s, t, u}, noorder];
    ML2 = createSQS[{v}, {w}, noorder];
    
    ML3 = createSQS[{p, q, r}, {t, u, w}, noorder];
    ML4 = createSQS[{p, q, r}, {s, t, w}, noorder];
    ML5 = createSQS[{p, q, r}, {s, u, w}, noorder];
    ML6 = createSQS[{p, q, r, v}, {s, t, u, w}, noorder];
    res = wick[ML1 ** ML2, {p, q, r, s, t, u, v, w}, wickopts];
    nor = 
    nor = Map[Distribute, nor];
    res === nor
    ,
    True
    ,
    TestID->"SeQuantTest_Multi_wick_3_1"
]*)