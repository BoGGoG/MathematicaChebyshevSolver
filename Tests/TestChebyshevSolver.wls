(* TEST CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript
(* ::Package:: *)

<<"../source/ChebyshevSolver.wl"

(*<<"../../MBI-Resurgence/source/DeriveEinsteinEquationsFollowingCartwright.wl"

<<"../../MBI-Resurgence/source/NumericalPart.wl"

<<"../../MathematicaPackagesGit/DEQSystemStrucutureVisualization/source/DEQSystemStructureVisualization.wl"*)

DEQ = 1 + f[x] + 3 g[x] f[x] f''[x] + 2 f''[x] f'[x] + 123 f'[x] == 0 ;
bc1 = f[0] == 0;
bc2 = f[1] == 0;
chebPoints = 25;
x0 = -3.;
x1 = 3.;

ChebyshevSolver`GetNthOrderTerm[DEQ, f, {x, 2}] == 3 g[x] f[x] f''[x] + 2 f''[x] f'[x]
