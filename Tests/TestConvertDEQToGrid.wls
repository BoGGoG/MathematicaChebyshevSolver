#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

(* Only analytic expressions *)

DEQ =  Sin[x] + x^2 f[x] + 2 f'[x] + f''[x] + Exp[Cos[x]] f'''[x]  == 0 ;
(* DEQ =  Sin[x] + x^2 f[x] + 2 f'[x] + f''[x]  == 0 ; *)
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];

DEQOnGrid = ConvertDEQToGrid[DEQ, f, {x, grid}];
DEQOnGridCompare = Map[GetCoefficientArray[DEQ, f, {x, #}, grid]&, {-1,0,1,2,3}];
Print[DEQOnGrid == DEQOnGridCompare];

(* including some functions as arrays on grid already *)

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

DEQ = 1 + x f[x] + 2 x Bs[x] f'[x] + f''[x] == 0;

Print[GetCoefficientArray[DEQ, f, {x, -1}, grid] == GetCoefficientArray[DEQ, f, {x, -1}, grid, {{Bs, BsArr}}]];
Print[GetCoefficientArray[DEQ, f, {x, 0}, grid] == GetCoefficientArray[DEQ, f, {x, 0}, grid, {{Bs, BsArr}}]];
Print[GetCoefficientArray[DEQ, f, {x, 2}, grid] == GetCoefficientArray[DEQ, f, {x, 2}, grid, {{Bs, BsArr}}]];
Print[GetCoefficientArray[DEQ, f, {x, 1}, grid, {{Bs, BsArr}}] == Table[BsArr[[i]] 2 grid[[i]], {i,1,Length@grid}]];
