(* USE CHEBYSHEV SOLVER *)
#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"
<<"../source/ChebyshevSolverTester.wl"

(* DEQOrder *)

<<"../source/ChebyshevSolver.wl"

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 2];

DEQ =  x^2 + x^2 f[x] + 2 f'[x] == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 1];

DEQ =  x^2 + x^2 f[x]  == 0 ;
Print[ChebyshevSolver`DEQOrder[DEQ, f, x] == 0];

(* Compare to NDSolve *)

<<"../source/ChebyshevSolver.wl"
<<"../source/ChebyshevSolverTester.wl"

DEQ =  Cos[x] x^2 - Sin[x] f[x] + 2 Exp[x] f'[x] == 0 ;
x0 = -2.;
x1 = 4 Pi;
bc1 = f[x0] == 0.;
chebPoints = 200;
plot = ChebyshevSolverTester`TestDEQ[{DEQ, bc1}, f, {x,x0,x1}, "GridPoints"->chebPoints];
Export["../TmpPlots/plot.png", plot];

DEQ =  x^2 + x^2 f[x] + 2 f'[x] + f''[x] == 0 ;
x0 = -1.;
x1 = 4 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == 1.;
chebPoints = 200;
plot = ChebyshevSolverTester`TestDEQ[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];
Export["../TmpPlots/plot.png", plot];
