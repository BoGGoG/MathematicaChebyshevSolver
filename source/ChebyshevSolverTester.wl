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

Options[GeneratePlots] = {"GridPoints" -> 50, "NumberOfDigits"->MachinePrecision};
GeneratePlots[DEQAndBCs__, f_, {x_, x0_, x1_}, OptionsPattern[]] := Block[
		{sol, grid, solInterpolated, solNDSolve, plotCheb, plotNDSolve, plotInterpolate},

	{sol, grid} = ChebyshevSolver`ChebyNDSolveRaw[DEQAndBCs, f, {x,x0,x1}, "GridPoints"->OptionValue["GridPoints"]];
	solInterpolated = ChebyshevSolver`ChebyNDSolve[DEQAndBCs, f, {x,x0,x1}, "GridPoints"->OptionValue["GridPoints"]];
	solNDSolve[y_] = f[y]/.NDSolve[DEQAndBCs, f, {x, x0, x1}][[1]];

	plotCheb = ListPlot[Thread[{grid, sol}], PlotLegends->LineLegend@{"ChebyNDSolveRaw"}];
	plotInterpolate = Plot[solInterpolated[x], {x,x0,x1},
		PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (interpolated)"}];
	plotNDSolve = Plot[solNDSolve[x], {x,x0,x1}, PlotStyle->{Red, Dashed}, PlotLegends->LineLegend@{"NDSolve"}];
	{plotCheb, plotInterpolate, plotNDSolve}
];

ShowPlots[{plotCheb_, plotNDSolve_, plotInterpolate_}, DEQAndBCs__, nGrid_, {x_, x0_, x1_}] := Block[
	{},
	Show[plotCheb, plotNDSolve, plotInterpolate,
		PlotLabel->"DEQ: "<>ToString[TraditionalForm@DEQAndBCs[[1]]]<>"\nBCs: "<>ToString[DEQAndBCs[[2;;]]]
		<>"\nGrid: ["<>ToString@x0<>","<>ToString@x1"]"<>", GridPoints: "<>ToString[nGrid]]
];


Options[TestDEQ] = {"GridPoints" -> 50, "NumberOfDigits"->MachinePrecision};
TestDEQ[DEQAndBCs__, f_, {x_, x0_, x1_}, OptionsPattern[]] := Block[
		{DEQ, bcs, y, plotCheb, plotNDSolve, plotInterpolate},

	PrintSetup[DEQAndBCs, {x, x0, x1}];
	{plotCheb, plotInterpolate, plotNDSolve} = GeneratePlots[DEQAndBCs, f, {x, x0, x1},
		"GridPoints"->OptionValue["GridPoints"]];
	ShowPlots[{plotCheb, plotNDSolve, plotInterpolate}, DEQAndBCs, OptionValue["GridPoints"], {x,x0,x1}]
];

CreateDirIfNotExists[path_] := Quiet[
  CreateDirectory[path]
, CreateDirectory::filex
];

Options[SavePlots] = {"Folder"->"../TmpPlots"};
SavePlots[plots_, OptionsPattern[]] := Block[{name, p, path},
	CreateDirIfNotExists[OptionValue["Folder"]];
	Print[Length@plots, " plots to export"];
	Do[
		name = "DEQ"<>ToString@p<>".png";
		path = FileNameJoin[{OptionValue["Folder"], name}];
		Print["Exporting plot ", path];
		Export[path, plots[[p]]];
	, {p, 1, Length@plots}];
];

TestDEQs[setups_] := Block[{plots},
	plots = Map[Apply[TestDEQ], setups];
	SavePlots[plots];
];

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["NumericalPart`*"], Head[#] === Symbol &]];
EndPackage[];

(* ToDo
	- check boundary conditions
*)
