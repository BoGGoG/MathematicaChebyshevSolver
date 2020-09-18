(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

(*<<"../../MBI-Resurgence/source/DeriveEinsteinEquationsFollowingCartwright.wl"

<<"../../MBI-Resurgence/source/NumericalPart.wl"

<<"../../MathematicaPackagesGit/DEQSystemStrucutureVisualization/source/DEQSystemStructureVisualization.wl"*)

DEQ = 1 + f[x] +  2 f''[x] == 0 ;
bc1 = f[0] == 0;
bc2 = f[1] == 0;
chebPoints = 25;
x0 = -3.;
x1 = 3.;

ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}]

(* Plot and Export for Viewing *)
plot = MatrixPlot[DerivMatrix];
Export["../TmpPlots/plot.png", plot]
