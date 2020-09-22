(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

(*<<"../../MBI-Resurgence/source/DeriveEinsteinEquationsFollowingCartwright.wl"
<<"../../MBI-Resurgence/source/NumericalPart.wl"
<<"../../MathematicaPackagesGit/DEQSystemStrucutureVisualization/source/DEQSystemStructureVisualization.wl"*)

<<"../source/ChebyshevSolver.wl"

DEQ = + x f[x] - 2 f'[x] + f''[x] == 0 ;
(* DEQ = f[x] + 2 f'[x] - f''[x] == 0 ; *)
(* DEQ = Sin[x]^2 f[x] + Exp[-x] f''[x] == 0 (* needs chebPoints = 100 *) *)
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 100;
{sol, grid} = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

sol2[y_] = f[y]/.NDSolve[{DEQ, bc1, bc2}, f, {x, x0, x1}][[1]];


(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{grid, sol}]];
plotNDSolve = Plot[sol2[x], {x,x0,x1}];
plot = Show[plotCheb, plotNDSolve];
Export["../TmpPlots/plot.png", plot]
(* Export["../TmpPlots/plot.png", DerivMatrix // MatrixForm] *)
