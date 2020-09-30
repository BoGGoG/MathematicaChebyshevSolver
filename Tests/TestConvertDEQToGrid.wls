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

(* including some functions as
