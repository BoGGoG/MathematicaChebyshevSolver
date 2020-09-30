#!/usr/bin/env wolframscript

(*
	ConvertDEQToGrid takes the DEQ and a grid and returns the coefficients
	as arrays as they are evaluated on the grid.
	A DEQ = s[x] + a[x] f[x] + b[x] f'[x] + ... gives {s[grid], a[grid], b[grid], ...}

	If one of the coefficient functions is problematic at zero, you can set
	"LimitPointIndex"->1 and at that point it will be treated more carefully.
	An example is x Log[x], which Mathematica cannot evaluate at zero, but it is actually
	finite there (0).

	ConvertDEQToGrid can also take parts of coefficients that are already on grid,
	just call it like
	ConvertDEQToGrid[DEQ, f, {x, grid}, OnGridFuncsAndValues], where
	OnGridFuncsAndValues = {{B, BArr}, ...}, where B is the function that is in the coefficients
	of DEQ and BArr are its values on the grid only.
	*)

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

DEQOnGrid = ConvertDEQToGrid[DEQ, f, {x, grid}, {{Bs, BsArr}}];

cn1Arr = GetCoefficientArray[DEQ, f, {x, -1}, grid];
c0Arr = GetCoefficientArray[DEQ, f, {x, 0}, grid];
c1Arr = GetCoefficientArray[DEQ, f, {x, 1}, grid, {{Bs, BsArr}}];
c2Arr = GetCoefficientArray[DEQ, f, {x, 2}, grid];
DEQOnGridCompare = {cn1Arr, c0Arr, c1Arr, c2Arr};
Print[DEQOnGrid == DEQOnGridCompare];

cn1ArrManually = Table[1, {i,1,Length@grid}];
c0ArrManually = Table[grid[[i]], {i,1,Length@grid}];
c1ArrManually = Table[2 grid[[i]] BsArr[[i]], {i,1,Length@grid}];
c2ArrManually = Table[1, {i,1,Length@grid}];
DEQOnGridManually = {cn1ArrManually, c0ArrManually, c1ArrManually, c2ArrManually};
Print[DEQOnGrid == DEQOnGridManually];
