#!/usr/bin/env wolframscript

<<"../source/ChebyshevSolver.wl"

DEQ =  Sin[x] + x^2 f[x] + 2 f'[x] + f''[x] + Exp[Cos[x]] f'''[x]  == 0 ;
x0 = 0.;
x1 = 2 Pi;
bc1 = f[x0] == 0.;
bc2 = f'[x0] == + 1.5;
chebPoints = 5;
{sol, {grid, deriv}} = ChebyshevSolver`ChebyNDSolveRaw[{DEQ, bc1, bc2}, f, {x,x0,x1}, "GridPoints"->chebPoints];

allTests = {};

(* GetCoefficients *)
coefficients = ChebyshevSolver`GetCoefficients[DEQ, f, {x,3}];
coefficients == {Sin[x], x^2, 2, 1, Exp[Cos[x]]}

(* coefficient arrays *)
coefficientsArrays = Map[GetCoefficientArray[DEQ, f, {x, #}, grid]&, Range[0,3]];
coeff1ArrCompare = coefficients[[1]]/.x->grid;
coeff2ArrCompare = coefficients[[2]]/.x->grid;
coeff3ArrCompare = ConstantArray[coefficients[[3]], Length@grid];
coeff4ArrCompare = ConstantArray[coefficients[[4]], Length@grid];
coeff5ArrCompare = coefficients[[5]]/.x->grid;

coeff2ArrCompare == coefficientsArrays[[1]]
coeff3ArrCompare == coefficientsArrays[[2]]
coeff4ArrCompare == coefficientsArrays[[3]]
coeff5ArrCompare == coefficientsArrays[[4]]

(* coefficient that is indep of f *)
indepCoeffArrCompare = (DEQ[[1]]/.f->(0&))/.x->grid;
indepCoeffArr = ChebyshevSolver`GetIndepCoeff[DEQ, f, x]/.x->grid;
indepCoeffArrCompare == indepCoeffArr

plot = ListPlot[Thread@{grid, arr}];

DEQOnGrid = ChebyshevSolver`ConvertDEQToGrid[DEQ, f, x, grid]

orderZeroOp = ChebyshevSolver`BuildDEQMatrixOrderNFromGridValues[DEQOnGrid[[2]], {grid, deriv}, 0]

DEQMatrixOp = ChebyshevSolver`BuildDEQMatrixFromGridValues[DEQOnGrid, {grid, deriv}]

DEQMatrixOp//MatrixForm


plot = ListPlot[Thread@{grid, DEQOnGrid[[2]]}];


Export["../TmpPlots/plot.png", plot] (* might have to create this folder before *)

Export["../TmpPlots/plot.png", DEQMatrixOp//MatrixForm] (* might have to create this folder before *)


DEQMatrixOp
