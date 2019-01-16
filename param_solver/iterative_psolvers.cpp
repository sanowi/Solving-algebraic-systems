//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: 
//%     Iterative methods for solving parametric interval linear systems
//%
//%     Contains:
//%     SpectralRadius, firstSolution, InnerIter
//%     CheckForZeros, Accurate, StrictlyContains, ResidIter
//%     ZVector, CMatrix
//%     Rump, ParamGaussSeidel, Krawczyk, MKolev, Jacobi
//%
//%##########################################################################
//%

#include <conio.h>
#include "../utils/stdafx.h"
#include "../interval/interval_base.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../param_solver/direct_solvers.h"
#include "../param_solver/hansenbliekrohn.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"

const real eps = 1.0e-3;
const real mi = spn(); //% smallest positive number

#define SIGN(x)		(((x) > 0) ? 1 : ((x) < 0) ? -1 : 0)
#define EPSX		1.0e-5 //% accuracy
#define MAXITER		500 //% maximal number of iterations

int
//%----------------------------------------------------------------------------
//% Computes the spectral radius of matrix B=Ac^{-1}*A, where A is a revised
//% affine matrix, and optionally prints B
//%----------------------------------------------------------------------------
SpectralRadius(int op, //% op = 0 - left pre-cond, op = 1 - right pre-cond
	int pr, //% print option
	const aafrmatrix& A, //% revised affine matrix
	dmatrix& B) //% radius of conditioned matrix A
{
	int n = A.num_rows();
	dmatrix R(n, n);
	aafrmatrix C(n, n), Rf(n, n);
	aafrvector b(n);

	mid(A, R); //% midpoint matrix
	if (pr == 1) {
		std::cout << "B=[";
	}
	if (inv(R)) { //% midpoint inverse
		if (op == 0) {
			C = R * A; //% left pre-conditioning
		}
		else {
			C = A * R; //% right pre-conditining
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = C(i, j).reduce().rad();
				if (pr == 1) {
					std::cout << B(i, j) << " ";
				}
			}
			if (pr == 1) {
				if (i < n - 1) {
					std::cout << ";" << std::endl;
				}
			}
		}
		if (pr == 1) {
			std::cout << "];" << std::endl << std::endl;
			std::cout << std::endl;
		}
		return Success;
	}
	return SingularMatrix;
}

real
//%----------------------------------------------------------------------------
//% Computes the spectral radius of matrix B=Ac^{-1}*A, where A is an
//% interval matrix, and optionally prints B
//%----------------------------------------------------------------------------
SpectralRadius(int op, //% op = 0 - left pre-cond, op = 1 - right pre-cond
	const imatrix& A) //% revised affine matrix
{
	int n = A.num_rows();
	dmatrix R(n, n), Cr(n, n);
	imatrix C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		if (op == 0) {
			C = R * A; //% left pre-conditioning
		}
		else {
			C = A * R; //% right pre-conditining
		}
		rad(C, Cr);
		return rhoSpectral(Cr);
	}
	return ISLIB_INFINITY;
} //% >>> spectralRadius <<<

void
//%----------------------------------------------------------------------------
//% >>> IterInner <<<
//% inner solution computed from a vector of reduced affine forms
//%----------------------------------------------------------------------------
IterInner(const aafrvector& x, //% p-solution (revised affine vector)
	dvector& innl, //% lower bound of IEH
	dvector& innu) //% upper bound of IEH
{
	int n = x.size();
	dvector	lam(n), ia(n);

	for (int i = 0; i < n; i++) {
		int ltemp = x[i].idx().size();
		lam[i] = 0.0;
		for (int k = 1; k < ltemp; k++) {
			lam[i] = lam[i] + ISLIB_ABS(x[i].cfs()[k]);
		}
		innl[i] = x[i].cfs()[0] - lam[i] + x[i].r();
		innu[i] = x[i].cfs()[0] + lam[i] - x[i].r();
	}
} //% >>> IterInner <<<

void
//%----------------------------------------------------------------------------
//% >>> IterInner <<<
//% inner solution computed from a matrix of revised affine forms
//% (for systems with multiple right-hand side)
//%---------------------------------------------------------------------------- 
IterInner(const aafrmatrix& x, //% p-solution (matrix of revised affine forms)
	dmatrix& innl, //% lower bound if inner solution
	dmatrix& innu) //% upper bound if inner solution
{
	int n = x.num_rows();
	int m = x.num_cols();
	dmatrix lam(n, m), ia(n, m);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int ltemp = x(i, j).idx().size();
			lam(i, j) = 0.0;
			for (int k = 1; k < ltemp; k++) {
				lam(i, j) = lam(i, j) + ISLIB_ABS(x(i, j).cfs()[k]);
			}
			innl(i, j) = x(i, j).cfs()[0] - lam(i, j) + x(i, j).r();
			innu(i, j) = x(i, j).cfs()[0] + lam(i, j) - x(i, j).r();
		}
	}
} //% >>> IterInner <<<

static void
//%----------------------------------------------------------------------------
//% >>> CheckForZeros <<<
//% The vectors x and y are successive approximations for the solution of 
//% a linear system of equations computed by iterative refinement. If 
//% a component of y is diminished by more than 'Factor', it is a good 
//% candidate for a zero entry. Thus, it is set to zero.
//%----------------------------------------------------------------------------
//% static function is "private" to its .c file (or its translation unit)
//%----------------------------------------------------------------------------
CheckForZeros(dvector& x, dvector& y)
{
	const real Factor = 1e+5;
	int n = x.size();
	assert(x.size() == y.size());

	for (int i = 0; i < n; i++) {
		if (ISLIB_ABS(y[i])*Factor < ISLIB_ABS(x[i])) y[i] = 0.0;
	}
} //% >>> CheckForZeros <<<

static bool
//%----------------------------------------------------------------------------
//% >>> Accurate <<<
//% The vectors x and y are successive iterates. The function returns true if
//% the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i
//% has about 12 correct decimals. If x_i or y_i vanishes, the relative error
//% is implicitly set to zero.
//%----------------------------------------------------------------------------
Accurate(dvector& x, dvector& y)
{
	const real Delta = 1e-12;   //% Relative error bound
	int i = 0, n = y.size();
	bool ok = false;
	real abs_yi;

	assert(x.size() == y.size());
	do {
		if (SIGN(x[i])*SIGN(y[i]) <= 0.0) //% Relative error set to 0
			ok = true;
		else {
			abs_yi = ISLIB_ABS(y[i]);
			ok = (ISLIB_ABS(y[i] - x[i]) <= Delta * abs_yi); //% Relative error > Delta
		}
		i++;
	} while (ok && (i < n));
	return ok;
} //% >>> Accurate <<<

void
//%----------------------------------------------------------------------------
//% >>> ResidIter <<<
//% Iteratively improves solution of the linear system Ax=b
//%----------------------------------------------------------------------------
ResidIter(const dmatrix& R, //% midpoint inverse
	const dmatrix& A, //% sysmte matrix
	const dvector& b, //% right-hand vector
	dvector& x) //% approximate improved solution
{
	int n = x.size(), k = 0;
	bool ok = false;
	dvector y(n), d(n);

	do {
		d = b - A * x; //% should be rounded probably
		y = x + R * d; //% should be rounded probably
		CheckForZeros(x, y);
		ok = Accurate(x, y);
		if (ok || k > MAXITER) return;
		x = y;
		k++;
	} while (true);
} //% >>> ResidIter <<<

static void
//%----------------------------------------------------------------------------
//% Computes R * M(p), where M(p) is omvector, i.e., vector of matrices
//%----------------------------------------------------------------------------
CMatrix(const dmatrix& R,
	const ivector& p,
	const omvector& M,
	imatrix& C)
{
	int n = M[0].num_cols(), K = p.size();
	omvector B(K + 1);

	for (int k = 0; k < K + 1; k++) {
		B[k] = R * M[k];
	}
	C = imatrix(n, n);
	C.fill_in(0.0);
	C = C + B[0];
	for (int k = 1; k < K + 1; k++) {
		C = C + B[k] * p[k - 1];
	}
}

static void
//%----------------------------------------------------------------------------
//% Computes R * (b(p) - M(p) * x0), where M(p) is omvector, i.e., 
//% vector of matrices and b(p) is vector of vectors
//%----------------------------------------------------------------------------
ZVector(const dvector& x0,
	const dmatrix& R,
	const ivector& p,
	const omvector& M,
	const ovector& b,
	ivector& z)
{
	int n = b[0].size(), K = p.size();

	z = ivector(n);
	z.fill_in(0.0);
	z = z + R * (b[0] - M[0] * x0);
	for (int k = 1; k < K + 1; k++) {
		z = z + R * (b[k] - M[k] * x0) * p[k - 1];
	}
}

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using revised affine forms.
//% Linear dependencies are passed through iterations.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% interval outer solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector z(n), y0(n), y1(n);
	aafrmatrix B(n, n), C(n, n);
	dvector bc(n), x0(n);
	dmatrix Ac(n, n), Bm(n, n), R(n, n);

	w = aafrvector(n);
	mid(A, Ac); //% midpoint matrix
	R = Ac; //% copy of midpoint matrix
	if (inv(R)) { //% inverse of midpoint matrix
		mid(b, bc); //% midpoint vector
		x0 = R * bc; //% midpoint solution
		ResidIter(R, Ac, bc, x0); //% residual correction
		C = R * A; //% preconditioning
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = (i == j) ? 1.0 - C(i, j) : -C(i, j);
			}
		}
		mag(reduce(B), Bm); //% |I-RA|
		if (rhoSpectral(Bm) >= 1.0) { return NotStronglyRegular; }
		z = R * (b - A * x0); //% residual correction and left-preconditioning
		y0 = z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (verification step)
		//%-------------------------------------------------------------
		while (true) {
			niter++;
			if (niter > MAXITER) { break; } //% maximal number of iterations exeeded
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; i++) {
				AAFR ykk = 0.0;
				for (int j = 0; j < n; j++) {
					ykk = ykk + B(i, j) * y1[j];
				}
				y1[i] = z[i] + ykk;
			}
			if (reduce(y0).set_strictly_contains(reduce(y1))) { 
				for (int i = 0; i < n; i++) {
					w[i] = x0[i] + y1[i];
				}
				return Success; //% verified solution has been found
			}
			y0 = y1;
		}
		if (dist1(y0, y1) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Rump <<<

int
//%---------------------------------------------------------------------------- 
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The method uses epresentation of parametric systems as
//% an affine linear combination of real matrices, where the parameters play
//% the role of the coefficients.
//% ---!!! Method used in global optimization !!!---
//%----------------------------------------------------------------------------
Rump(const ivector& p, //% vector of interval parameters
	const omvector& oA, //% matrix of left hand dependencies
	const ovector& ob, //% vector of right hand dependencies
	ivector& w, //% interval solution vector
	int& stepCount) //% number of iterations
{
	if (!oA[0].is_square()) return NonSquareMatrix;

	int	n = oA[0].num_rows(), K = p.size();
	dvector	bm(n), x0(n), yylb(n), yyub(n), pm(K), y0r(n);
	dmatrix	Cm(n, n), R(n, n), I(n, n);
	imatrix	A(n, n), B(n, n), C(n, n);
	ivector	b(n), z(n), y0(n), y1(n), yk(n), yy(n), imi(n);

	R = mid(oA, p); //% midpoint matrix
	bm = mid(ob, p); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		x0 = R * bm; //% midpoint solution
		CMatrix(R, p, oA, C); //% C = R*A(p)
		Unit(I);
		B = I - C;
		mig(B, Cm); // Cm = <C>
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		ZVector(x0, R, p, oA, ob, z); //% z = R*(b(p)-A(p)*x0)
		y0 = z;
		stepCount = 0;
		//% Main iteration loop
		while (true) {
			stepCount++;
			if (stepCount > MAXITER) { return Failed; }
			blow(y0, eps, mi);
			y1 = y0;

			for (int i = 0; i < n; i++) {
				yk = z + (B * y1);
				y1[i] = yk[i];
			}
			if (y0.set_strictly_contains(yk)) {
				w = x0 + yk;
				return Success;
			}
			y0 = yk;
		}
	}
	return SingularMatrix;
} //% >>> Rump with oA and ob <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using revised affine forms.
//% Linear dependencies are passed through iterations.
//% !!!Right!!! preconditioning
//%----------------------------------------------------------------------------
RumpII(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% interval outer solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector z(n), y0(n), y1(n);
	aafrmatrix B(n, n), C(n, n);
	dvector bc(n), x0(n);
	dmatrix Ac(n, n), Bm(n, n), R(n, n);

	w = aafrvector(n);
	mid(A, Ac); //% midpoint matrix
	R = Ac; //% copy of midpoint matrix
	if (inv(R)) { //% inverse of midpoint matrix
		mid(b, bc); //% midpoint vector
		x0 = R * bc; //% midpoint solution
		ResidIter(R, Ac, bc, x0); //% residual correction
		C = A * R; //% preconditioning
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = (i == j) ? 1.0 - C(i, j) : -C(i, j);
			}
		}
		mag(reduce(B), Bm); //% |I-RA|
		if (rhoSpectral(Bm) >= 1.0) { return NotStronglyRegular; }
		z = b - A * x0; //% residual correction and left-preconditioning
		y0 = z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (verification step)
		//%-------------------------------------------------------------
		while (true) {
			niter++;
			if (niter > MAXITER) { break; } //% maximal number of iterations exeeded
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; i++) {
				AAFR ykk = 0.0;
				for (int j = 0; j < n; j++) {
					ykk = ykk + B(i, j) * y1[j];
				}
				y1[i] = z[i] + ykk;
			}
			if (reduce(y0).set_strictly_contains(reduce(y1))) {
				y1 = R * y1;
				for (int i = 0; i < n; i++) {
					w[i] = x0[i] + y1[i];
				}
				return Success; //% verified solution has been found
			}
			y0 = y1;
		}
		if (dist1(y0, y1) < EPSX) {
			y1 = R * y1;
			for (int i = 0; i < n; i++) {
				w[i] = x0[i] + y1[i];
			}
			return NonVerifiedSolution;
		}
		y1 = R * y1;
		for (int i = 0; i < n; i++) {
			w[i] = x0[i] + y1[i];
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Rump <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% for systems with *** Multiple Right-hand Side (MHRS) ***
//% Self-verified fixed point iteration for solving parametric interval 
//% systems with multiple righ-hand side. The method is based on Krawczyk 
//% iteration.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& b, //% revised affine right-hand side matrix
	aafrmatrix& w, //% revised affine solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows(), m = b.num_cols();
	dmatrix	Bc(n, m), X0(n, m);
	dmatrix	Cm(n, n), R(n, n), I(n, n);
	aafrmatrix	B(n, n), C(n, n);
	aafrmatrix	Z(n, m), Y0(n, m), Y1(n, m), yk(n, m), yy(n, m), D(n, m);
	aafrmatrix Ax(n, m);
	AAFR ykk;

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% inverse of the midpoint matrix
		mid(b, Bc); //% midpoint vector
		X0 = R * Bc; //% midpoint solution
		C = R * A; //% left preconditioning
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - C(i, j);
			}
		}
		mag(reduce(B), Cm);
		if (rhoSpectral(Cm) >= 1.0) { return NotStronglyRegular; }
		Z = R * (b - A * X0); //% right preconditioning and residual correction
		Y0 = Z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		while (true) {
			niter++;
			if (niter > MAXITER) { break; }
			blow(Y0, eps, mi); //% blowing up is used to obtain verified solution
			Y1 = Y0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					AAFR ykk = 0.0;
					for (int k = 0; k < n; k++) {
						ykk = ykk + B(i, k) * Y1(k, j);
					}
					Y1(i, j) = Z(i, j) + ykk;
				}
			}
			if (StrictlyContains(reduce(Y0), reduce(Y1))) { //% solution found
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						w(i, j) = X0(i, j) + Y1(i, j);
					}
				}
				return Success; //% verified solution has been found
			}
			Y0 = Y1;
		}
		if (dist1(Y0, Y1) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Rump for MHRS <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% self verified method with pre-conditioning using the midpoint inverse
//% The iterative scheme of the method is as follows:
//% x=D^{-1}(b-(L+U)x)
//% It origines from the following decompsition of the system matrix:
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
//%---------------------------------------------------------------------------- 
ParamGaussSeidel(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows();
	aafrvector w1(n), c(n);
	aafrmatrix C(n, n);
	dmatrix R(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		c = R * b; //% O(n^2*K)
		//w.fill_in(0.0);
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			niter++;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < i; j++) { //% lower triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; j++) { //% upper triangular
					if (C(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
						s = s + C(i, j) * w1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				w[i] = (c[i] - s) / C(i, i); //% O(K*log(K))
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(w1, w) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> ParamGaussSeidel <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% self verified method with right pre-conditioning with midpoint inverse
//%----------------------------------------------------------------------------
ParamGaussSeidelRPrec(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows();
	aafrvector w1(n);
	aafrmatrix C(n, n);
	dmatrix	R(n, n), Am(n, n);

	mid(A, R);
	if (inv(R)) {
		C = A * R; //% right pre-conditioning
		w.fill_in(0.0);
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			niter++;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			for (int i = 0; i < n; i++) { // all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
				AAFR s = 0.0;
				for (int j = 0; j < i; j++) { //% lower triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; j++) { //% upper triangular
					if (C(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
						s = s + C(i, j) * w1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				w[i] = (b[i] - s) / C(i, i); //% O(K*log(K))
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				w = R * w; //% final solution
				return Success;  //% verified solution
			}
		} while (true);
		if (distH(reduce(w1), reduce(w)) < EPSX) {
			w = R * w; //% final non-verified solution
			return NonVerifiedSolution;
		}
		w = R * w;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> ParamGaussSeidelRPrec <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% with double (left- and right) pre-conditioning using the midpoint inverse
//%----------------------------------------------------------------------------
ParamGaussSeidelIII(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	const dmatrix& L, //% left preconditioner
	const dmatrix R, //% right prconditioner
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows(), err;
	aafrvector w1(n), c(n);
	aafrmatrix C(n, n);
	dmatrix	R2(n, n), Am(n, n);
	ivector xp(n);
	
	//%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//% PRZEMYSLEC!!!
	err = ParametricDirectMethodPIII2(A, b, L, R, xp);
	//% initial enclosure of the right pre-conditioned system
	if (err != Success) return err;
	//% transform initial interval enclosure into revised affine form
	//% with midpoint of interval as central value and radius of interval
	//% as radius of accumulative error, no noise variable is introduced
	w.fill_in(0.0);
	for (int i = 0; i < n; i++) { // O(n)
		w[i] = w[i] + xp[i];
	}
	C = L * A * R; //% right pre-conditioning
	c = L * b;
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	err = Success;
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) {
			break;
		}
		w1 = w;
		for (int i = 0; i < n; i++) { //% all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) { //% lower triangular
				if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
					s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) { //% upper triangular
				if (C(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
					s = s + C(i, j) * w1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			w[i] = (c[i] - s) / C(i, i); //% O(K*log(K))
		}
	} while (dist1(w, w1) > EPSX);
	w = R * w; //% final solution
	return err;
} //% >>> ParamGaussSeidelII <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel <<<
//% iteration with pre-conditioning for solving systems with 
//% multiple right-hand side
//%----------------------------------------------------------------------------
ParamGaussSeidel(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& b, //% revised right-hand side vector
	aafrmatrix& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows(), m = b.num_cols();
	aafrmatrix w1(n, m), c(n, m);
	aafrmatrix C(n, n);
	dmatrix R(n, n);

	w = aafrmatrix(n, m);
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		c = R * b; //% O(n^2*K)
		w.fill_in(0.0);
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			niter++;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < m; k++) {
					AAFR s = 0.0;
					for (int j = 0; j < i; j++) { //% lower triangular
						if (C(i, j).reduce() != 0.0 && w(j, k).reduce() != 0.0) {
							s = s + C(i, j) * w(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					for (int j = i + 1; j < n; j++) { //% upper triangular
						if (C(i, j).reduce() != 0.0 && w1(j, k).reduce() != 0.0) {
							s = s + C(i, j) * w1(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					w(i, k) = (c(i, k) - s) / C(i, i); //% O(K*log(K))
				}
			}
			if (StrictlyContains(reduce(w1), reduce(w))) {
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(w1, w) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> ParamGaussSeidel PROBA <<<

int
//%----------------------------------------------------------------------------
//% >>> Krawczyk iteration <<< 
//% without residual correction; the starting vector must be computed by the
//% method without residual correction 
//%----------------------------------------------------------------------------
Krawczyk(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector v(n), y0(n), y1(n);
	aafrmatrix V(n, n), B(n, n);
	dmatrix	R(n, n), Br(n, n);
	w = aafrvector(n);
	mid(A, R);
	if (inv(R)) {
		V = R * A;
		v = R * b;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = (i == j) ? 1.0 - V(i, j) : -V(i, j);
			}
		}
		mag(reduce(B), Br);
		real rs = rhoSpectral(Br);
		if (rs >= 1.0) { return NotStronglyRegular; }
		//% transform initial interval enclosure into revised affine form
		//% with midpoint of interval as central value and radius of interval
		//% as radius of accumulative error, no noise variable is introduced
		//y0.fill_in(0.0);
		y0 = v;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) {
				break;
			}
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < n; j++) {
					if (B(i, j).reduce() != 0.0 && y0[j].reduce() != 0.0) {
						s = s + B(i, j) * y0[j];
					}
				}
				y0[i] = v[i] + s; //% update entry
			}
			if (reduce(y1).set_strictly_contains(reduce(y0))) {
				w = y0;
				return Success; //% verified solution has been found
			}
		} while (true);
		if (distH(y0, y1) < EPSX) {
			w = y0;
			return NonVerifiedSolution;
		}
		w = y0;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Krawczyk <<<

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration wihtout residual correction
//% the starting vector must be computed by the
//% method without residual correction 
//%----------------------------------------------------------------------------
KrawczykII(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows(), err;
	aafrvector v(n), y0(n), y1(n);
	aafrmatrix C(n, n), B(n, n);
	dmatrix	R(n, n), I(n, n), Br(n, n);
	ivector	w0(n);

	w = aafrvector(n);
	//% initial solution
	err = HansenBliekRohnII0(A, b, w0, R);
	if (err != Success) return err;
	C = A * R; //% right preconditioning of A
	Unit(I);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B(i, j) = I(i, j) - C(i, j); //% (I-V(e))
		}
	}
	mag(reduce(B), Br);
	real rs = rhoSpectral(Br);
	if (rs >= 1.0) {
		return NotStronglyRegular;
	}
	//% transform initial interval enclosure into revised affine form
	//% with midpoint of interval as central value and radius of interval
	//% as radius of accumulative error, no noise variable is introduced
	y0.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		y0[i] = y0[i] + w0[i]; //% starting point
	}
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) { return MaxIterExceeded; }
		y1 = y0;
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				if (B(i, j).reduce() != 0.0 && y0[j].reduce() != 0.0) {
					s = s + B(i, j) * y0[j];
				}
			}
			y0[i] = b[i] + s; //% update entry
		}
	} while (distH(y0, y1) > EPSX);
	w = R * y0;
	return Success;
} //% >>> KrawczykII <<<

int
//%----------------------------------------------------------------------------
//% >>> Krawczyk iteration <<<
//% described in the paper with Milan. The method is adapted to solve systems
//% with multiple right-hand side.
//%----------------------------------------------------------------------------
Krawczyk(const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& bb, //% revised right-hand side vector
	aafrmatrix& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = bb.num_cols();
	real eps = 0.01, mi = 1.0e-15, rs;
	aafrmatrix C(n, n), B(n, n), c(n, m), y0(n, m), y1(n, m);
	dmatrix R(n, n), Br(n, n);
	
	w = aafrmatrix(n, m);
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		//% left pre-conditioning
		C = R * A;
		c = R * bb;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = (i == j) ? 1.0 - C(i, j) : -C(i, j); //% matrix (I-V(e))
			}
		}
		mag(reduce(B), Br);
		rs = rhoSpectral(Br);
		if (rs >= 1.0) return NotStronglyRegular;
		//% transform initial interval enclosure into revised affine form
		//% with midpoint of interval as central value and radius of interval
		//% as radius of accumulative error, no noise variable is introduced
		y0.fill_in(0.0);
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; i++) { //% all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
				for (int k = 0; k < m; k++) {
					AAFR s = 0.0;
					for (int j = 0; j < n; j++) {
						if (B(i, j).reduce() != 0.0 && y0(j, k).reduce() != 0.0) {
							s = s + B(i, j) * y0(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					y0(i, k) = c(i, k) + s; //% O(K*log(K))
				}
			}
			if (StrictlyContains(reduce(y1), reduce(y0))) {
				w = y0;
				return Success; //% verified solution has been found
			}

		} while (true);
		if (dist1(y0, y1) < EPSX) {
			w = y0;
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Krawczyk <<<

int
//%----------------------------------------------------------------------------
//% >>> Modified Kolev's Iteration method <<< 
//% without residual correction. The method is based on revised affine forms
//%----------------------------------------------------------------------------
MKolev(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector c(n), v(n), v1(n);
	aafrmatrix C(n, n), C0(n, n);
	dmatrix A0(n, n), R(n, n);

	mid(A, R);
	if (inv(R)) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C0(i, j) = A(i, j) - A(i, j).cfs()[0];
			}
		}
		//% left pre-conditioning
		C = R * C0;
		c = R * b;
		//% transform initial interval enclosure into revised affine form
		//% with midpoint of interval as central value and radius of interval
		//% as radius of accumulative error, no noise variable is introduced
		v.fill_in(0.0);
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) { break; }
			blow(v, eps, mi); //% blowing up is used to obtain verified solution
			v1 = v;
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < n; j++) {
					if (C(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + C(i, j) * v[j];
					}
				}
				v[i] = c[i] - s; //% update entry
			}
			if (reduce(v1).set_strictly_contains(reduce(v))) {
				w = v;
				return Success; //% verified solution has been found
			}
		} while (true); 
		if (dist1(v1, v) < EPSX) {
			w = v;
			return NonVerifiedSolution;
		}
		w = v;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> MKolev <<<

int
//%----------------------------------------------------------------------------
//% Modified version of Kolev's method
//% computing with reduced affine forms
//% with !!!right!!! pre-conditioning
//%----------------------------------------------------------------------------
MKolevII(const aafrmatrix& A, //% revied affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w, //% p-solution (revised affine matrix)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector d(n), x0f(n), v(n), v1(n);
	aafrmatrix C(n, n), B(n, n), Af(n, n);
	cvector	xpidx(1);
	dvector x0(n), xpcfs(1);
	dmatrix R(n, n);
	ivector xp(n);

	mid(A, R);
	if (inv(R)) {
		B = A * R; //% right pre-conditioning of system matrix
		d = b;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C(i, j) = B(i, j) - B(i, j).cfs()[0];
			}
		}
		v = b;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			blow(v, eps, mi);
			v1 = v;
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < n; j++) {
					if (C(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + C(i, j) * v[j];
					}
				}
				v[i] = d[i] - s;
			}
			if (reduce(v1).set_strictly_contains(reduce(v))) {
				w = R * v;
				return Success;
			}
		} while (true);
		if (dist1(v, v1) < EPSX) {
			w = R * v; //% final solution
			return NonVerifiedSolution;
		}
		w = R * v;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> MKolevII <<<

int
//%----------------------------------------------------------------------------
//% Modified Kolev's Iteration method 
//% with residual correction and right pre-conditioning
//% based on revised affine forms
//%----------------------------------------------------------------------------
MKolevRCII(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations 
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector d(n), v(n), v1(n);
	dvector bm(n), x0(n);
	ivector xp(n);
	aafrmatrix C(n, n), C0(n, n);
	dmatrix	A0(n, n), R(n, n);

	mid(b, bm); //% midpoint vector
	err = ParametricDirectMethodPII(A, b, xp, R); //% xp - initial enclosure (right pre-cond)
	if (err != Success) return err;
	x0 = R * bm; //% midpoint solution
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C0(i, j) = A(i, j) - A(i, j).cfs()[0];
		}
	}
	C = C0 * R; //% right pre-conditioning
	d = (b - A * x0); //% residual correction
					  //% transform initial interval enclosure into revised affine form
	v.fill_in(0.0);
	for (int i = 0; i < n; i++) { //% O(n)
		v[i] = v[i] + xp[i];
	}
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		v1 = v;
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) { //% use already updated entries
				if (C(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
					s = s + C(i, j) * v[j];
				}
			}
			v[i] = d[i] - s;
		}
	} while (dist1(v1, v) > EPSX);
	v = R * v;
	for (int i = 0; i < n; i++) {
		w[i] = (x0[i] + v[i]);
	}
	return Success;
} //% >>> MKolevRCII <<<

int
//%----------------------------------------------------------------------------
//% >>> Modified Kolev's Iteration method <<<
//% with residual correction
//% based on revised affine forms
//% this is almost the same as Krawczyk iteration
//%---------------------------------------------------------------------------- 
MKolevRC(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter)	//% out: number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector d(n), v(n), v1(n);
	aafrmatrix C(n, n), C0(n, n);
	dvector bm(n), x0(n);
	dmatrix A0(n, n), R(n, n);

	mid(b, bm); //% midpoint vector mid(b)
	mid(A, R);
	if (inv(R)) {
		x0 = R * bm; //% midpoint solution
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C0(i, j) = A(i, j) - A(i, j).cfs()[0];
			}
		}
		d = R * (b - A * x0); //% residual correction
		C = R * C0;
		//% transform initial interval enclosure into revised affine form
		v.fill_in(0.0);
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			v1 = v;
			blow(v, eps, mi);
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < n; j++) { //% use already updated entries
					if (C(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + C(i, j) * v[j];
					}
				}
				v[i] = d[i] - s;
			}
			if (reduce(v1).set_strictly_contains(reduce(v))) {
				for (int i = 0; i < n; i++) {
					w[i] = (x0[i] + v[i]);
				}
				return Success; //% verified solution has been found
			}

		} while (true);
		if (dist1(v1, v) < EPSX) {
			for (int i = 0; i < n; i++) {
				w[i] = (x0[i] + v[i]);
			}
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> MKolevRC <<<

int
//%----------------------------------------------------------------------------
//% >>> Jacobi iteration <<<
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%---------------------------------------------------------------------------- 
Jacobi(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& x0,	//% p-solution (revised affine vector)
	int& niter)	//% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector x1(n), c(n);
	aafrmatrix C(n, n);
	dmatrix R(n, n);

	mid(A, R);
	if (inv(R)) {
		x0.fill_in(0.0);
		C = R * A;
		c = R * b;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			x1 = x0;
			blow(x0, eps, mi);
			for (int i = 0; i < n; i++) {
				AAFR s = 0.0;
				for (int j = 0; j < n; j++) {
					if (i != j) {
						if (C(i, j).reduce() != 0.0 && x1[j].reduce() != 0.0) {
							s = s + C(i, j) * x1[j];
						}
					}
				}
				x0[i] = (c[i] - s) / C(i, i);
			}
			if (reduce(x1).set_strictly_contains(reduce(x0))) {
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(x0, x1) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Jacobi <<<

int
//%----------------------------------------------------------------------------
//% Parametric Gauss-Seidel iteration
//% Iterative scheme: x=D^{-1}(b-(L+U)x)
//% Origin: Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% the difference is that 
//%----------------------------------------------------------------------------
ParamGaussSeidel0(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand vector
	aafrvector& x0,	//% interval solution
	int& niter)	//% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector bb(n), x1(n);
	cvector xpidx(1);
	dvector xinf(n), xsup(n), xpcfs(1);
	ivector	xp(n), xi(n), xii(n);
	aafrmatrix AA(n, n);
	dmatrix R(n, n);

	//% initial enclosure obtained using ISM
	err = ParametricDirectMethod(A, b, xp, R); // O(n^3*K)
	if (err != Success) return err;
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		x0[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//% left pre-conditioning with the midpoint inverse
	AA = R * A; //% O(n^3*K)
	bb = R * b; //% O(n^2*K)
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		xi = reduce(x1); //% O(n)
		for (int i = 0; i < n; i++) { //% all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) {
				if (AA(i, j).reduce() != 0.0 && x0[j].reduce() != 0.0) {
					s = s + AA(i, j) * x0[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) {
				if (AA(i, j).reduce() != 0.0 && x1[j].reduce() != 0.0) {
					s = s + AA(i, j) * x1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			xii[i] = ((bb[i] - s) / AA(i, i)).reduce();
			xii[i] = interval(ISLIB_MAX(xi[i].inf(), xii[i].inf()), ISLIB_MIN(xi[i].sup(), xii[i].sup()));
			xpcfs[0] = xii[i].mid();
			x0[i] = AAFR(xpidx, xpcfs, xii[i].rad());
		}
	} while (distH(xi, xii) > EPSX);
	return Success;
} //% >>> ParamGaussSeidel0 <<<

int
//%----------------------------------------------------------------------------
//% Succesive overrelaxation method
//% Iterative method (preconditioning with minpoint inverse)
//%---------------------------------------------------------------------------- 
SuccessiveOverrelaxation(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand vector
	aafrvector& x0, //% interval solution
	real wgt, //% relaxation parameter
	int& niter)	//% number of iterations
{
	int n = A.num_rows(), err;
	cvector xpidx(1);
	dvector xpcfs(1);
	ivector xp(n);
	dmatrix R(n, n);
	aafrvector x1(n), bb(n);
	aafrmatrix AA(n, n), Rf(n, n);

	err = ParametricDirectMethod(A, b, xp, R); //% initial enclosure
	if (err != Success) return err;
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		x0[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	AA = R * A;
	bb = R * b;
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		for (int i = 0; i < n; i++) {
			AAFR s1 = 0.0;
			AAFR s2 = 0.0;
			for (int j = 0; j < i; j++) {
				s1 = s1 + AA(i, j) * x0[j];
			}
			for (int j = i + 1; j < n; j++) {
				s2 = s2 + AA(i, j) * x1[j];
			}
			x0[i] = (1.0 - wgt)*x1[i] + wgt * (bb[i] - s1 - s2) / AA(i, i);
		}
	} while (dist1(x0, x1) > EPSX);
	return Success;
} //% >>> SuccessiveOverrelaxation <<<