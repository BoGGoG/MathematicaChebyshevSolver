(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"
<<"../source/ChebyshevSolverTester.wl"

(* Test DEQOrder *)
(* DEQ order prints the order of the ordinary differential equation *)

<<"../source/ChebyshevSolver.wl"

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 2];

DEQ =  x^2 + x^2 f[x] + 2 f'[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 1];

DEQ =  x^2 + x^2 f[x]  == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 0];

<<"../source/ChebyshevSolver.wl"
<<"../source/ChebyshevSolverTester.wl"

(* Test GetCoefficientArray *)
DEQ = s[x] + a[x] f[x] + b[x] f'[x] + c[x] f''[x] == 0;
x0 = 0.;
x1 = 2 Pi;
chebPoints = 5;
{grid, deriv} = ChebyshevSolver`ChebyshevSetup[chebPoints, "Intervall"->{x0,x1}];
Print[ ChebyshevSolver`GetCoefficientArray[DEQ, f, {x, -1}, grid] == Map[s, grid] ];
Print[ ChebyshevSolver`GetCoefficientArray[DEQ, f, {x, 0}, grid] == Map[a, grid] ];
Print[ ChebyshevSolver`GetCoefficientArray[DEQ, f, {x, 1}, grid] == Map[b, grid]] ;
Print[ ChebyshevSolver`GetCoefficientArray[DEQ, f, {x, 2}, grid] == Map[c, grid]] ;


(* Test
	- ChebyNDSolveRaw: Returns grid and solution on grid
	- ChebyNDSolveRaw: Returns interpolated solution (as function)
*)

setupList = {};

DEQ =  f[x] + f''[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 50;
setup = {{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints};
AppendTo[setupList, setup];

DEQ =  Cos[x] x^2 - Sin[x] f[x] + 2 Exp[x] f'[x] == 0 ;
x0 = -2.;
x1 = 4 Pi;
bc1 = f[x0] == 0.;
chebPoints = 50;
setup = {{DEQ, bc1}, f, {x,x0,x1}, "GridPoints"->chebPoints};
AppendTo[setupList, setup];

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
x0 = -1.;
x1 = 4 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 100;
setup = {{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints};
AppendTo[setupList, setup];

DEQ = f[x] + 2 f'[x] + f''[x] == 0 ;
x0 = -1.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 100;
setup = {{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints};
AppendTo[setupList, setup];

DEQ = f[x] + Exp[x] f''[x] == 0 ;
x0 = -1.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 100;
setup = {{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints};
AppendTo[setupList, setup];

DEQ = x Log[x] + x Log[x] f[x] + f'[x] + f''[x] == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 100;
setup = {{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints, "LimitPointIndex"->1};
AppendTo[setupList, setup];

TestDEQs[setupList]
