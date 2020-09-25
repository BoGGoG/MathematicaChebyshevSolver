(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
(* DEQ = f[x] + Exp[-x] f''[x] == 0; *)
(* DEQ = f[x] + 2 f'[x] - f''[x] == 0 ; *)
(* DEQ = Sin[x]^2 f[x] + Exp[-x] f''[x] == 0 (* needs chebPoints = 100 *) *)
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 50;
{sol, {grid, deriv}} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

{solInterpolated, dsolInterpolated} = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];
sol2[y_] = f[y]/.NDSolve[{DEQ, bc1, bc2}, f, {x, x0, x1}][[1]];
(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{grid, sol}], PlotLegends->LineLegend@{"ChebyNDSolveRaw"}];
plotNDSolve = Plot[sol2[x], {x,x0,x1}, PlotStyle->{Red, Dashed}, PlotLegends->LineLegend@{"NDSolve"}];
plotInterpolate = Plot[solInterpolated[x], {x,x0,x1}, PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (interpolated)"}];

plot = Show[plotCheb, plotNDSolve, plotInterpolate];
Export["../TmpPlots/plot.png", plot] (* might have to create this folder before *)

dplotInterpolate = Plot[solInterpolated'[x], {x,x0,x1}, PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (interpolated)"}];
dplotMatrixMult = Plot[dsolInterpolated[x], {x,x0,x1}, PlotStyle->{Opacity[0.4,Blue]}, PlotLegends->LineLegend@{"ChebyNDSolve (matrix mult)"}];
dplot = Show[dplotInterpolate, dplotMatrixMult];
Export["../TmpPlots/dplot.png", dplot] (* might have to create this folder before *)

{1,2}, {3,4}, {5,6}

{{{1}, 3, 5}, {{2}, 4, 6}}
