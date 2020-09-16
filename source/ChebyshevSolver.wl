BeginPackage["ChebyshevSolver`"];

Unprotect["ChebyshevSolver`*"];
ClearAll["ChebyshevSolver`*", "ChebyshevSolver`Private`*"];

Begin["`Private`"];

examplefunction[]:=Print["Hello World!"]

StripArgument[expr_, arg_] := expr /. f_[xx___, arg, yy___] -> f[xx, yy];

CGLGrid[x0_, L_, n_Integer /; n > 1] :=
 x0 + 1/2 L (1 -  Cos[\[Pi] Range[0, n - 1]/(n - 1)])

Scan[SetAttributes[#, {Protected, ReadProtected}]&,
     Select[Symbol /@ Names["NumericalPart`*"], Head[#] === Symbol &]];

End[];
EndPackage[];
