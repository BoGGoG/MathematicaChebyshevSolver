#!/usr/bin/env wolframscript
(* ::Package:: *)

<<"../source/ChebyshevSolver.wl"

(*<< JavaGraphics`*)

(* Importing Differential Equations *)
DEQ = f[x] + f''[x] ;
chebPoints = 25;
x0 = -3.;
x1 = 3.;

(* DEQ1 is for Ss, so find coefficients of Ss *)

fieldAndDerivatives = {f[x], f'[x], f''[x]};
coeffs = Coefficient[DEQ, fieldAndDerivatives];

(* Chebyshev Setup *)
{cGrid, DerivMatrix} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0, x1}];

(* Build Operator *)
DEQMatrixOperator = coeffs[[1]] IdentityMatrix[chebPoints+1] + coeffs[[2]] DerivMatrix + coeffs[[3]] DerivMatrix.DerivMatrix;

bc1 = Table[If[i==1,1,0], {i,1,chebPoints+1}];
bc2 = Table[If[i==chebPoints+1,1,0], {i,1,chebPoints+1}];
rhs = Table[0, chebPoints+1];
rhs[[1]] = 1;
rhs[[-1]] = 2;
DEQMatrixOperator[[1]] = bc1;
DEQMatrixOperator[[-1]] = bc2;

sol = LinearSolve[DEQMatrixOperator, rhs];

(*sol[x_] = Total@Map[(ChebyshevT[#-1, x] eigVec[[1]][[#]] &), Range[1,chebPoints]]*)

sol2[y_] = f[y]/.NDSolve[{DEQ==0, f[x0] == 1, f[x1] == 2}, f, {x, x0, x1}][[1]]

(* Plot and Export for Viewing *)
plotCheb = ListPlot[Thread[{cGrid, sol}]];
plotNDSolve = Plot[sol2[x], {x,x0,x1}];
plot = Show[plotCheb, plotNDSolve];
Export["../TmpPlots/plot.png", plot]
