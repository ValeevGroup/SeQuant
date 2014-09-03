(* Mathematica Test File *)

(* Test Function for wick function in MultiConfiguration *)

(*

normal contraction, using all occupied orbitals

*)


(* wick between one and one operator *)
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
	deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] SQM[OHead["\[Lambda]", indexSymm[-1]],
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
 deltaIndex[{particleIndex["q", 
     particleSpace[occupied]]}, {particleIndex["r", 
     particleSpace[occupied]]}] SQS[
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