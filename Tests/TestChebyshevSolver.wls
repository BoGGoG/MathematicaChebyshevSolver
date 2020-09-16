#!/usr/bin/env wolframscript
(* ::Package:: *)

<<"../source/ChebyshevSolver.wl"

<<"../../MBI-Resurgence/source/DeriveEinsteinEquationsFollowingCartwright.wl"

<<"../../MBI-Resurgence/source/NumericalPart.wl"

<<"../../MathematicaPackagesGit/DEQSystemStrucutureVisualization/source/DEQSystemStructureVisualization.wl"

<< JavaGraphics`

(* Importing Differential Equations *)
DEQs = Import[FileNameJoin[{"..","Equations", "EFEsFinal.m"}]];
fieldsLambdas = {Bs, Ss, dplusSs, dplusBs, As};
fields = Through[fieldsLambdas[t,u]];
fieldsStriped = ChebyshevSolver`Private`StripArgument[fields, t];
fieldsAllDerivatives = Map[AllUpToNthDerivatives[#, 2]&, fieldsStriped];

maxDerivative = 2;

ShowDEQSystemStructure[DEQs, fields]

(* DEQ1 is for Ss, so find coefficients of Ss *)

SsFieldAndDerivs = {Ss[t,u], D[Ss[t,u],u], D[Ss[t,u],u,u]};

SsCoefficients = Coefficient[DEQs[[1]], SsFieldAndDerivs];

(* Chebyshev Setup *)
{cGrid, DerivMatrix} = ChebyshevSolver`Private`ChebyshevSetup[0., 1., 25];





(* Plot and Export for Viewing *)
plot = MatrixPlot[DerivMatrix];
Export["../TmpPlots/plot.png", plot]
