(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"
<<"../source/ChebyshevSolverTester.wl"

DEQ =  x^2 + x^2 f[x] + 2 Exp[x] f'[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
chebPoints = 50;
(* Tests *)
plot = ChebyshevSolverTester`TestDEQ[{DEQ, bc1}, f, {x,x0,x1}, "GridPoints"->chebPoints]
Export["../TmpPlots/plot.png", plot]

(* DEQ SOLVING *)

(* First order DEQs *)
DEQ =  x^2 + x^2 f[x] + 2 Exp[x] f'[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
chebPoints = 50;

{sol, grid} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1}, f, {x,x0,x1}, "GridPoints"->chebPoints];
solInterpolated = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1}, f, {x,x0,x1}, "GridPoints"->chebPoints];

solNDSolve[y_] = f[y]/.NDSolve[{DEQ, bc1}, f, {x, x0, x1}][[1]];

(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{grid, sol}], PlotLegends->LineLegend@{"ChebyNDSolveRaw"}];
plotNDSolve = Plot[solNDSolve[x], {x,x0,x1}, PlotStyle->{Red, Dashed}, PlotLegends->LineLegend@{"NDSolve"}];
plotInterpolate = Plot[solInterpolated[x], {x,x0,x1}, PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (interpolated)"}];
plot = Show[plotCheb, plotNDSolve, plotInterpolate];
Export["../TmpPlots/plot.png", plot]

(* DEQOrder *)

<<"../source/ChebyshevSolver.wl"

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 2];

DEQ =  x^2 + x^2 f[x] + 2 f'[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 1];

DEQ =  x^2 + x^2 f[x]  == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 0];
