(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"
x0 = 0.;
x1 = 1.;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];

(* Test 1 *)

coeff = x Sin[x] + x Bs[x] + Ss[x];
coeffFunc = x Sin[x] + x (x Exp[-x] ) + Cos[x];
SsArr = Cos[x] /. x->grid;
BsArr = x Exp[-x] /. x->grid;
gridVals = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}, {Ss, SsArr}}];
gridValsCompare = coeffFunc /. x->grid

gridVals == gridValsCompare

(* Test 2 *)


SsFunc[x_] := Cos[x];
BsFunc[x_] := x Exp[-x];
AsFunc[x_] := Exp[x];
SsArr = Cos[x] /. x->grid;
BsArr = x Exp[-x] /. x->grid;
AsArr = Exp[x] /. x->grid;
coeff = x Sin[x] + Sin[x] Bs[x] + Ss[x] + x As[x];
coeffFunc = x Sin[x] + Sin[x] BsFunc[x] + SsFunc[x] + x AsFunc[x];

gridVals = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}, {Ss, SsArr}, {As, AsArr}}]
gridValsCompare = coeffFunc /. x->grid

gridVals == gridValsCompare


(* Test 3 *)


SsFunc[x_] := Cos[x];
BsFunc[x_] := x Exp[-x];
AsFunc[x_] := Exp[x];
SsArr = Cos[x] /. x->grid;
BsArr = x Exp[-x] /. x->grid;
AsArr = Exp[x] /. x->grid;
coeff = x Sin[x] + Ss[x] Bs[x] +  x As[x];
coeffFunc = x Sin[x] + SsFunc[x] BsFunc[x] + x AsFunc[x];

gridVals = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}, {Ss, SsArr}, {As, AsArr}}]
gridValsCompare = coeffFunc /. x->grid

gridVals == gridValsCompare

(* What about points where one cannot plug in the point directly? *)
coeff = Cos[x] + x Log[x] + Bs[x];
gridVals = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}}, "LimitPointIndex"->1];
gridValsIndet = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}}]
specialPoint = coeff /. Bs[x]->BsArr[[1]] /. x Log[x] -> 0 /. x->grid[[1]]

Print[gridVals[[2;;]] == gridValsIndet[[2;;]]];
Print[gridVals[[1]] == specialPoint];
