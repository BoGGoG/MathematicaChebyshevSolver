(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

(*<<"../../MBI-Resurgence/source/DeriveEinsteinEquationsFollowingCartwright.wl"
<<"../../MBI-Resurgence/source/NumericalPart.wl"
<<"../../MathematicaPackagesGit/DEQSystemStrucutureVisualization/source/DEQSystemStructureVisualization.wl"*)

<<"../source/ChebyshevSolver.wl"

DEQ = f[x] +  f''[x] == 0 ;
bc1 = f[0.] == 0;
bc2 = f'[0.] == 1;
chebPoints = 15;
x0 = 0.;
x1 = 3.14.;
{DEQMatrixOperator, rhs, grid} = ChebyshevSolver`ChebyNDSolve[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

sol = LinearSolve[DEQMatrixOperator, rhs];
sol2[y_] = f[y]/.NDSolve[{DEQ==0, f[x0] == 1, f[x1] == 2}, f, {x, x0, x1}][[1]]

(* Plot and Export for Viewing *)
plot = ListPlot[Thread@{grid, sol}];
Export["../TmpPlots/plot.png", plot]
(* Export["../TmpPlots/plot.png", DerivMatrix // MatrixForm] *)
