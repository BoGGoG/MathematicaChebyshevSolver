(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

x0 = 0.;
x1 = 1.;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];


Bs[x_] := 1 + Exp[-x] Sin[x];
BsArr = Bs[x] /. x->grid;

coeffFunctions = x Sin[x] + x Bs[x];

coeff = x Sin[x] + x Bss[x];

<<"../source/ChebyshevSolver.wl"
ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{"Bs", "BsArr"}, {"Ss", "SsArr"}}]

coeff2 = coeff /.Bss[x]->Bss


Table[coeff2 /. {x->grid[[i]], Bss -> BsArr[[i]]}, {i, 1, Length@grid}]

Table[coeffFunctions /. x->grid]


(* ---------------------- *)

<<"../source/ChebyshevSolver.wl"
x0 = 0.;
x1 = 1.;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];
coeff = x Sin[x] + x Bs[x] + Ss[x];
coeffFunc = x Sin[x] + x (x Exp[-x] ) + Cos[x];
SsArr = Cos[x] /. x->grid;
BsArr = x Exp[-x] /. x->grid;
ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}, {Ss, SsArr}}]

coeffFunc /. x->grid

Transpose[{{Bs, BsArr}, {Ss, SsArr}}][[1]]


(* ---------------------- *)
(* Special Evaluate on Grid *)
<<"../source/ChebyshevSolver.wl"

x0 = 0.;
x1 = 1.;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];
coeff = x Log[x] + Cos[x] + Bs[x] + Ss[x];
BsArr = x Exp[-x] /. x->grid;
SsArr = Cos[x] /. x->grid;
gridVals = ChebyshevSolver`EvaluateOnGrid[coeff, {x, grid}, {{Bs, BsArr}, {Ss, SsArr}}]

SpecialEvaluateOnGrid[coeff, {x, grid}, 1, {{Bs, BsArr}, {Ss, SsArr}}]

coeff /. Bs[x]->BsArr[[1]] /. Ss[x]->SsArr[[1]] /. Log[x]->0 /. x->grid[[1]]

ChebyshevSolver`SpecialEvaluate[coeff, {x, 0.}, {{Bs, BsArr}, {Ss, SsArr}}]
