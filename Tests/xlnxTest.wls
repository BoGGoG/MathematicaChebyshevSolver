#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

DEQ =   x Log[x] f[x] + (2*x)*f'[x] + x Exp[x] f''[x]  == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 50;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];
DEQOnGrid = ChebyshevSolver`ConvertDEQToGrid[DEQ, f, x, grid, "LimitPointIndex"->1];

solND[x_] = f[x]/.NDSolve[{DEQ, f[0] == 0, f'[0] == 1.5}, f, {x,x0, x1}][[1]];


{sol0, {grid, deriv}} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1, bc2}, f, {x,x0,x1},
	"GridPoints"->chebPoints, "LimitPointIndex"->1];
{sol, {grid, deriv}} = ChebyshevSolver`ChebyRawInputSolve[{DEQOnGrid, {bc1, bc2}}, {grid, deriv}];


plotND = Plot[solND[x], {x,x0,x1}];
plot0 = ListPlot[Thread@{grid, sol0}, PlotMarkers->"OpenMarkers", PlotStyle->Blue];
plot = ListPlot[Thread@{grid, sol}, PlotStyle->Red];
Export["../TmpPlots/plot2.png", Show[plot0, plot, plotND]] (* might have to create this folder before *)
