#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

DEQ =  Sin[x] + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 50;
{sol, {grid, deriv}} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

dfCoeff= ChebyshevSolver`GetNthOrderCoeff[DEQ, f, {x, 1}]
ddfCoeff= ChebyshevSolver`GetNthOrderCoeff[DEQ, f, {x, 2}]
fCoeff= ChebyshevSolver`GetNthOrderCoeff[DEQ, f, {x, 0}]
IndepCoeff= ChebyshevSolver`GetNthOrderCoeff[DEQ, f, {x, -1}]
