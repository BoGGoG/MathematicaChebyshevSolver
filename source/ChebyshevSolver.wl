BeginPackage["ChebyshevSolver`"];

Unprotect["ChebyshevSolver`*"];
ClearAll["ChebyshevSolver`*", "ChebyshevSolver`Private`*"];

(* Begin["`Private`"]; *)

StripArgument[expr_, arg_] := expr /. f_[xx___, arg, yy___] -> f[xx, yy];

(* build chebyshev grid *)
ChebyshevPoints[nz_, OptionsPattern["NumberOfDigits"->MachinePrecision]] := Block[{numberOfDigits},
	numberOfDigits = OptionValue["NumberOfDigits"];
	N[Table[Cos[(r \[Pi])/nz],{r,0,nz}], numberOfDigits]
];

ChebDerivMatrixDiag[nz_, chebyshevPoints_] :=
	Block[{firstDiag, lastDiag},
	firstDiag = (2(nz+1-1)^2+1)/6;
	lastDiag = -firstDiag;
	DiagonalMatrix[Join[{firstDiag}, Map[-#/(2(1-#^2))&, chebyshevPoints[[2;;-2]]], {lastDiag}]]
];

ChebDerivMatrixOffDiag[nz_, chebyshevPoints_] := Table[
	If[i!=j,(-1)^(i+j)/(chebyshevPoints[[i]]-chebyshevPoints[[j]]) If[i==1||i==nz+1,2,1]/If[j==1||j==nz+1,2,1],0],
	{i,nz+1},{j,nz+1}];

ChebyshevPointsIfNotGiven[nz_,chebPointsGiven_, numberOfDigits_:MachinePrecision] :=  If[TrueQ[chebPointsGiven], ChebyshevPoints[nz, "NumberOfDigits"->numberOfDigits], chebPointsGiven];

(* build chebyshev grid derivative matrix *)
ChebDerivMatrix[nz_, OptionsPattern[{"ChebyshevPoints"->True,"NumberOfDigits"->MachinePrecision}]] := Block[{chebyshevPoints,firstDiag, lastDiag, diagMatrix,offDiag},

	chebyshevPoints = ChebyshevPointsIfNotGiven[nz, OptionValue["ChebyshevPoints"], OptionValue["NumberOfDigits"]];

	diagMatrix = ChebDerivMatrixDiag[nz, chebyshevPoints];
	offDiag = ChebDerivMatrixOffDiag[nz, chebyshevPoints];

	diagMatrix + offDiag
];

TransformChebIntervall[chebyshevPoints_, {a_,b_}] := chebyshevPoints*(a-b)/2+(a+b)/2;
(* using fact that it's symmetric to it's mid *)

TransformDCheb[DCheb_, {a_, b_}] := DCheb*2/(a-b);

(* build cheby grid and derivative matrix
returns {grid, matrix} *)
ChebyshevSetup[nz_, OptionsPattern[{"NumberOfDigits"->MachinePrecision, "Intervall"->{1,-1}}]]:=Block[
		{numberOfDigits, a,b, chebyshevPoints, DCheb},
	numberOfDigits = OptionValue["NumberOfDigits"];
	{a,b} = OptionValue["Intervall"];chebyshevPoints = ChebyshevPoints[nz, "NumberOfDigits"->numberOfDigits];DCheb =  ChebDerivMatrix[nz,"ChebyshevPoints"-> chebyshevPoints];
	If[{a,b}!={1,-1},
	chebyshevPoints = TransformChebIntervall[chebyshevPoints, {a,b}];
	DCheb = TransformDCheb[DCheb, {a,b}]];

	{chebyshevPoints, DCheb}
];

(* given cheby derivative matrices and coefficients of the ODE, build operator L
	such that L f = c f *)
BuildDEQMatrixOperator[coeffs_, deriv_] := Block[{nGrid, derivTerms, n},
	nGrid = Length@deriv;
	derivTerms = Sum[coeffs[[n]] MatrixPower[deriv,n-1], {n,2,Length@coeffs}];
	DEQMatrixOperator = coeffs[[1]] IdentityMatrix[nGrid] + derivTerms
]

HasDerivQ[expr_] := Not@FreeQ[expr, Derivative];

(* ApplyDirichletBC: not an easy task:
	- f[p] == y:
		- find grid point closest to p
		- set row in op corresponding to this grid point to be {0,...,1,0,0,...}
		- give a "source" term {0,...,y,0,0,...}
		- then op == source has f[p]==y as the row where p is
	- f'[p] == y:
*)
ApplyDirichletBC[DEQOperator_, bc_, {x_, x0_, x1_}] := Block[{operator, pos, f, p, nGrid},
	{pos, bcVal} = bc /. f_[p_] == y_ -> {p, y};
	nGrid = Length@DEQOperator;
	operator = DEQOperator;
	bcRow = Table[0, nGrid];
	rhs = Table[0, nGrid];

	If[Not@MemberQ[{x0,x1}, pos], Print["ERROR in ApplyDirichletBC, BC "<>ToString[bc]<>" is not at x0 or x1, this is not implemented yet!"]];

	If[pos == x0,
		rhs[[1]] = bcVal;
		bcRow[[1]] = 1;
		operator[[1]] = bcRow;
	];
	If[pos == x1,
		rhs[[-1]] = bcVal;
		bcRow[[-1]] = 1;
		operator[[-1]] = bcRow;
	];
	{operator, rhs}
]

(* assuming two boundary conditions, we only work with second order DEQs *)
AddBoundaryCond[DEQOperator_, boundaryConditions_, {x_, x0_, x1_}] := Block[{operator, rhs1, rhs2},
	operator = DEQOperator;
	If[Not@HasDerivQ[boundaryConditions[[1]]],
		{operator, rhs1} = ApplyDirichletBC[operator, boundaryConditions[[1]], {x,x0,x1}];,
		{operator, rhs1} = ApplyNeumannBC[operator, boundaryConditions[[1]], {x,x0,x1}];];

	If[Not@HasDerivQ[boundaryConditions[[2]]],
		{operator, rhs2} = ApplyDirichletBC[operator, boundaryConditions[[2]], {x,x0,x1}];,
		{operator, rhs2} = ApplyNeumannBC[operator, boundaryConditions[[2]], {x,x0,x1}];];

	{operator, rhs1+rhs2}
];

Options[ChebyNDSolve] = {"GridPoints" -> 25, "NumberOfDigits"->MachinePrecision};

(* only for up to second order ordinary DEQ *)
ChebyNDSolve[DEQAndBCs__, f_, {x_,x0_,x1_}, OptionsPattern[]] := Block[
		{DEQ, BCs, nGrid, funcAndDerivs, coeffs, DEQMatrixOperator, rhs},
	DEQ = DEQAndBCs[[1]][[1]];
	BCs = DEQAndBCs[[2;;]];
	nGrid = OptionValue["GridPoints"];

	{grid, deriv} = ChebyshevSetup[nGrid, "Intervall"->{x0,x1}, "NumberOfDigits"->OptionValue["NumberOfDigits"]];
	funcAndDerivs = Map[Derivative[#][f][x]&, {0,1,2}];
	coeffs = Coefficient[DEQ, funcAndDerivs];

	DEQMatrixOperator = BuildDEQMatrixOperator[coeffs, deriv];
	{DEQMatrixOperator, rhs} = AddBoundaryCond[DEQMatrixOperator, BCs, {x,x0,x1}];
	Print[rhs];
	DEQMatrixOperator
];

GetNthOrderTerm[DEQ_, f_, {x_, n_}] := Select[DEQ[[1]], Not[FreeQ[#, Derivative[n][f][x]]] &];


(* END OF FUNCTIONS *)

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["NumericalPart`*"], Head[#] === Symbol &]];

(*End[];*)
EndPackage[];
