(* Mathematica Test File *)


(* normal ordering for excitation operator with rank 1 *)
Test [
	operator = mSQS[normalOrder[False], 
 		SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  		particleIndex["q", particleSpace[occupied], indexType[ann]]]];
  	normalOrderedForm[operator]
  	,
  	SQM[OHead["\[Lambda]", indexSymm[-1]], 
  particleIndex["q", particleSpace[occupied], indexType[bra]], 
  particleIndex["p", particleSpace[occupied], indexType[ket]]] + 
 SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  particleIndex["q", particleSpace[occupied], indexType[ann]]]
  ,
  TestID->"SeQuantTest_Multi_normalOrder_1"
]


(* normal ordering for excitation operator with rank 2 *)
Test [
	operator = mSQS[normalOrder[False], 
 SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  particleIndex["r", particleSpace[occupied], indexType[cre]], 
  particleIndex["s", particleSpace[occupied], indexType[ann]], 
  particleIndex["q", particleSpace[occupied], indexType[ann]]]];
  
  normalOrderedForm[operator]
	,
	-SQM[OHead["\[Lambda]", indexSymm[-1]], 
    particleIndex["q", particleSpace[occupied], indexType[bra]], 
    particleIndex["r", particleSpace[occupied], indexType[ket]]] SQM[
   OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] + 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] SQM[
   OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]] + 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
  particleIndex["q", particleSpace[occupied], indexType[bra]], 
  particleIndex["s", particleSpace[occupied], indexType[bra]], 
  particleIndex["p", particleSpace[occupied], indexType[ket]], 
  particleIndex["r", particleSpace[occupied], indexType[ket]]] + 
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["s", particleSpace[occupied], indexType[bra]], 
   particleIndex["r", particleSpace[occupied], indexType[ket]]] SQS[
   particleIndex["p", particleSpace[occupied], indexType[cre]], 
   particleIndex["q", particleSpace[occupied], indexType[ann]]] - 
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
 SQM[OHead["\[Lambda]", indexSymm[-1]], 
   particleIndex["q", particleSpace[occupied], indexType[bra]], 
   particleIndex["p", particleSpace[occupied], indexType[ket]]] SQS[
   particleIndex["r", particleSpace[occupied], indexType[cre]], 
   particleIndex["s", particleSpace[occupied], indexType[ann]]] + 
 SQS[particleIndex["p", particleSpace[occupied], indexType[cre]], 
  particleIndex["r", particleSpace[occupied], indexType[cre]], 
  particleIndex["s", particleSpace[occupied], indexType[ann]], 
  particleIndex["q", particleSpace[occupied], indexType[ann]]]
  ,
  TestID->"SeQuantTest_Multi_normalOrder_2"
]