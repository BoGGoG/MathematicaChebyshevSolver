BeginPackage["ChebyshevSolverTester`", "ChebyshevSolver`"];

Unprotect["ChebyshevSolverTester`*"];
ClearAll["ChebyshevSolverTester`*", "ChebyshevSolverTester`Private`*"];

TestDEQ[{DEQAndBCs__}, f_, {x_, x0_, x1_}, OptionsPattern[]] := Block[{},
	0
];








Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["NumericalPart`*"], Head[#] === Symbol &]];
EndPackage[];
