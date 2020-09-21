(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

(* TEST 1 *)
DEQ = f[x] + 2 f'[x] +  f''[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 15;

{sol, grid} = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];
sol2[y_] = f[y]/.DSolve[{DEQ, bc1, bc2}, f, {x, x0, x1}][[1]]


(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{grid, sol}]];
plotNDSolve = Plot[sol2[x], {x,x0,x1}];
plot = Show[plotCheb, plotNDSolve];
Export["../TmpPlots/plot.png", plot]


(* TEST 2 *)

DEQ = f[x] + 2 f'[x] + Exp[x] f''[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 15;

{sol, grid} = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

sol2[y_] = f[y]/.NDSolve[{DEQ, bc1, bc2}, f, {x, x0, x1}][[1]]


(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{grid, sol}]];
plotNDSolve = Plot[sol2[x], {x,x0,x1}];
plot = Show[plotCheb, plotNDSolve];
Export["../TmpPlots/plot.png", plot]
