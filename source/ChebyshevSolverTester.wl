BeginPackage["ChebyshevSolverTester`", "ChebyshevSolver`"];

Unprotect["ChebyshevSolverTester`*"];
ClearAll["ChebyshevSolverTester`*", "ChebyshevSolverTester`Private`*"];

PrintSetup[DEQAndBCs__, {x_, x0_, x1_}] := Block[{DEQ, bcs},
	DEQ = DEQAndBCs[[1]];
	bcs = DEQAndBCs[[2;;]];
	Print["DEQ: ", DEQ];
	Print["Intervall: [", x0, ",", x1, "]"];
	Print["boundary conditions: ", bcs];
];

Options[TestDEQ] = {"GridPoints" -> 50, "NumberOfDigits"->MachinePrecision};
TestDEQ[DEQAndBCs__, f_, {x_, x0_, x1_}, OptionsPattern[]] := Block[
		{DEQ, bcs, y, plotCheb, plotNDSolve, plotInterpolate},
	PrintSetup[DEQAndBCs, {x, x0, x1}];

	{sol, grid} = ChebyshevSolver`ChebyNDSolveRaw[DEQAndBCs, f, {x,x0,x1}, "GridPoints"->OptionValue["GridPoints"]];
	solInterpolated = ChebyshevSolver`ChebyNDSolve[DEQAndBCs, f, {x,x0,x1}, "GridPoints"->OptionValue["GridPoints"]];
	solNDSolve[y_] = f[y]/.NDSolve[DEQAndBCs, f, {x, x0, x1}][[1]];

	plotCheb = ListPlot[Thread[{grid, sol}], PlotLegends->LineLegend@{"ChebyNDSolveRaw"}];
	plotNDSolve = Plot[solNDSolve[x], {x,x0,x1}, PlotStyle->{Red, Dashed}, PlotLegends->LineLegend@{"NDSolve"}];
	plotInterpolate = Plot[solInterpolated[x], {x,x0,x1},
		PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (interpolated)"}];
	Show[plotCheb, plotNDSolve, plotInterpolate,
		PlotLabel->"DEQ: "<>ToString[TraditionalForm@DEQAndBCs[[1]]]<>"\nBCs: "<>ToString[DEQAndBCs[[2;;]]]
		<>"\nGrid: ["<>ToString@x0<>","<>ToString@x1"]"<>", GridPoints: "<>ToString[OptionValue["GridPoints"]]]
];








Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["NumericalPart`*"], Head[#] === Symbol &]];
EndPackage[];
