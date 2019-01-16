#pragma once

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems based
//% on monotonicity approach
//%----------------------------------------------------------------------------
int MonotPsolution(int (*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), const aafrmatrix& A, const aafrvector& b, ivector& ww);
int MonotPsolutionMRHS(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), int(*ff)(const aafrmatrix&, const aafrmatrix&, aafrmatrix&, int&), const aafrmatrix& A, const aafrvector& b, ivector& ww);
int MonotStandardI(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), const aafrmatrix& A, const aafrvector& b, ivector& ww);
int MonotAugmentedD(int(*f)(const aafrmatrix&, const aafrvector&, ivector&), const aafrmatrix& A, const aafrvector& b, ivector& ww);
int MonotAugmentedI(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), const aafrmatrix& A, const aafrvector& b, ivector& ww);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems based
//% on monotonicity approach. Methods use representation of a matrix
//% as an affine combination of matrices. Not guaranteed results.
//%----------------------------------------------------------------------------
void Monot(const ivector& p, const omvector& M, const ovector& b, ivector& ww);
void MonotPsolution(int ii, const ivector& p, const omvector& M, const ovector& b, ivector& pmin, ivector& pmax);
void MonotPLP(int minmax, const dvector& c, const ivector& p, const omvector& M, const ovector& b, interval& ww);