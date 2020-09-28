#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

DEQ =  x Log[x] + x Log[x] f[x] + x Log[x] f'[x] + x Log[x] f''[x]  == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 50;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];
DEQOnGrid = ChebyshevSolver`ConvertDEQToGrid[DEQ, f, x, grid, "LimitPointIndex"->1];

{sol0, {grid, deriv}} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints,
	"LimitPointIndex"->1];

{sol, {grid, deriv}} = ChebyshevSolver`ChebyRawInputSolve[{DEQOnGrid, {bc1, bc2}}, {grid, deriv}];

plot0 = ListPlot[Thread@{grid, sol0}, PlotMarkers->"OpenMarkers", PlotStyle->Blue];
plot = ListPlot[Thread@{grid, sol}, PlotStyle->Red];
Export["../TmpPlots/plot2.png", Show[plot0, plot]] (* might have to create this folder before *)

sol0

ClearAll[func]


func[f_, {x0_, x1_}] := Block[{fnew,y, x},
	fnew[y_] := Switch[y,
		x0, Limit[f[x], x->x0],
		_, f[y]];
	 Map[fnew, {x0, x1}]
]

f := (# Log[#])&;

func[f, {0,1}]
