(* Mathematica Test File *)

(* Test Function for Contraction function in MultiConfiguration *)

(*

normal contraction, using all occupied orbitals

*)


(* contraction between one and one operator *)
Test[
	SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
	wickopts = {
   		fullContract -> False,
   		noCoincidences -> False,
   		doSums -> True,
   		doReindex -> True
   	};
	ML1 = SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
 particleIndex["q", particleSpace[occupied], indexType[ann]]];
	ML2 = SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
 particleIndex["s", particleSpace[occupied], indexType[ann]]];
	contractSQS[ ML1, ML2, wickopts] 
	,
	CR[1, SQS[particleIndex["p", particleSpace[occupied], indexType[cre]],
    particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[ann]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["s", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["s", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 1] +
  CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]], 1]
	,
	TestID->"SeQuantTest_Multi_contractSQS_1_1"
]

(* contraction between one and two operator *)
Test[
	SeQuantVacuum = SeQuantVacuumChoices["MultiConfiguration"];
	wickopts = {
   		fullContract -> False,
   		noCoincidences -> False,
   		doSums -> True,
   		doReindex -> True
   	};
	ML1 = SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
 particleIndex["q", particleSpace[occupied], indexType[ann]]];
	ML2 = SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
 particleIndex["s", particleSpace[occupied], indexType[cre]], 
 particleIndex["u", particleSpace[occupied], indexType[ann]], 
 particleIndex["t", particleSpace[occupied], indexType[ann]]];
	contractSQS[ ML1, ML2, wickopts] 
	,
	CR[1, SQS[particleIndex["p", particleSpace[occupied], indexType[cre]],
    particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Eta]", indexSymm[0]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Eta]", indexSymm[0]], 
     particleIndex["q", particleSpace[occupied], indexType[bra]], 
     particleIndex["s", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["u", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Eta]", indexSymm[0]], 
     particleIndex["q", particleSpace[occupied], indexType[bra]], 
     particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["t", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["t", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]], 
   particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["u", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 1] +
  CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["u", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]], 
   particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["t", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
     particleIndex["t", particleSpace[occupied], indexType[bra]], 
     particleIndex["p", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 1] +
  CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["t", particleSpace[occupied], indexType[bra]], 
   particleIndex["u", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["s", particleSpace[occupied], indexType[cre]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[-SQM[OHead["\[Eta]", indexSymm[0]], 
     particleIndex["q", particleSpace[occupied], indexType[bra]], 
     particleIndex["s", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]], 1] +
  CR[-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 
  SQS[particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]]] + 
 CR[SQM[OHead["\[Eta]", indexSymm[0]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
    OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["t", particleSpace[occupied], indexType[bra]], 
    particleIndex["u", particleSpace[occupied], indexType[bra]], 
    particleIndex["p", particleSpace[occupied], indexType[ket]], 
    particleIndex["s", particleSpace[occupied], indexType[ket]]], 1] +
  CR[SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["t", particleSpace[occupied], indexType[bra]], 
   particleIndex["u", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]], 
   particleIndex["s", particleSpace[occupied], indexType[ket]]], 1]
	,
	TestID->"SeQuantTest_Multi_contractSQS_1_2"
]
