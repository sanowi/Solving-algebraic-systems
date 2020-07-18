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
#include "../utils/svd.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../preconditioning/precond_spectral.h"

const real eps = 1.0e-3;
const real mi = spn(); //% smallest positive number

#define SIGN(x)		(((x) > 0) ? 1 : ((x) < 0) ? -1 : 0)
#define EPSX		1.0e-5 //% accuracy
#define MAXITER		500 //% maximal number of iterations

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
	dvector	lam{ n };

	for (int i = 0; i < n; ++i) {
		int ltemp = x[i].idx().size();
		lam[i] = 0.0;
		for (int k = 1; k < ltemp; ++k) {
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
//% for systems with multiple right-hand side
//%---------------------------------------------------------------------------- 
IterInner(const aafrmatrix& x, //% p-solution (matrix of revised affine forms)
	dmatrix& innl, //% lower bound if inner solution
	dmatrix& innu) //% upper bound if inner solution
{
	int n = x.num_rows(), m = x.num_cols();
	dmatrix lam{ n, m };

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			int ltemp = x(i, j).idx().size();
			lam(i, j) = 0.0;
			for (int k = 1; k < ltemp; ++k) {
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

	for (int i = 0; i < n; ++i) {
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
	const real delta = 1e-12;   //% Relative error bound
	int i = 0, n = y.size();
	bool ok = false;
	real abs_yi;

	assert(x.size() == y.size());
	do {
		if (SIGN(x[i])*SIGN(y[i]) <= 0.0) //% Relative error set to 0
			ok = true;
		else {
			abs_yi = ISLIB_ABS(y[i]);
			ok = (ISLIB_ABS(y[i] - x[i]) <= delta * abs_yi); //% Relative error > Delta
		}
		++i;
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
	dvector y{ n }, d{ n };

	do {
		d = b - A * x; //% should be rounded probably
		y = x + R * d; //% should be rounded probably
		CheckForZeros(x, y);
		ok = Accurate(x, y);
		if (ok || k > MAXITER) return;
		x = y;
		++k;
	} while (true);
} //% >>> ResidIter <<<

static void
//%----------------------------------------------------------------------------
//% Computes R * M(p), where M(p) is omvector, i.e., vector of matrices
//%----------------------------------------------------------------------------
cpmatrix(const dmatrix& R,
	const ivector& p,
	const omvector& M,
	imatrix& C)
{
	int n = M[0].num_cols(), K = p.size();
	omvector B{ K + 1 };
	C = imatrix{ n, n };
	C.fill_in(0.0);
	for (int k = 0; k < K + 1; ++k) {
		B[k] = R * M[k];
	}
	C = C + B[0];
	for (int k = 1; k < K + 1; ++k) {
		C = C + B[k] * p[k - 1];
	}
}//% >>> CMatrix <<<

static void
//%----------------------------------------------------------------------------
//% Computes R * (b(p) - M(p) * x0), where M(p) is omvector, i.e., 
//% vector of matrices and b(p) is vector of vectors
//%----------------------------------------------------------------------------
zvector(const dvector& xc, //% midpoint solution
	const dmatrix& R, //% inverse midpoint
	const ivector& p, //% vector of parameters
	const omvector& M, //% system matrux
	const ovector& b, //% system vector
	ivector& z)
{
	int n = b[0].size(), K = p.size();

	z = ivector{ n };
	z.fill_in(0.0);
	z = z + R * (b[0] - M[0] * xc);
	for (int k = 1; k < K + 1; ++k) {
		z = z + R * (b[k] - M[k] * xc) * p[k - 1];
	}
} //% >>> ZVector <<<

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
	real &rho,
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector z{ n }, y0{ n }, y1{ n };
	aafrmatrix B{ n, n }, C{ n, n };
	dvector bc{ n }, xc{ n };
	dmatrix Ac{ n, n }, Cm{ n, n }, R{ n, n };

	w = aafrvector{ n };
	mid(A, Ac); //% midpoint matrix
	R = Ac; //% copy of midpoint matrix
	if (inv(R)) { //% inverse of midpoint matrix
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		ResidIter(R, Ac, bc, xc); //% residual correction
		B = R * A; //% left preconditioning
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				C(i, j) = (i == j) ? 1.0 - B(i, j) : -B(i, j);
			}
		}
		mag(reduce(C), Cm); //% |I-RA|
		if ((rho = rhoSpectral(Cm)) >= 1.0) { return NotStronglyRegular; }
		z = R * (b - A * xc); //% residual correction and left-preconditioning
		y0 = z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (verification step)
		//%-------------------------------------------------------------
		while (true) {
			++niter;
			if (niter > MAXITER) { break; } //% maximal number of iterations exeeded
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) {
					s = s + C(i, j) * y0[j];
				}
				y0[i] = z[i] + s;
			}
			if (reduce(y1).set_strictly_contains(reduce(y0))) {
				w = xc + y0;
				return Success; //% verified solution has been found
			}
		}
		if (dist1(y0, y1) < EPSX) {
			w = xc + y0;
			return NonVerifiedSolution;
		}
		w = xc + y0;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Rump <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using revised affine forms.
//% Linear dependencies are passed through iterations.
//% !!!Right!!! preconditioning
//%----------------------------------------------------------------------------
RumpRCond(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% interval outer solution
	real &rho,
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector z{ n }, y0{ n }, y1{ n };
	aafrmatrix B{ n, n }, C{ n, n };
	dvector bc{ n }, xc{ n };
	dmatrix Ac{ n, n }, Bm{ n, n }, R{ n, n };

	w = aafrvector{ n };
	mid(A, Ac); //% midpoint matrix
	R = Ac; //% copy of midpoint matrix
	if (inv(R)) { //% inverse of midpoint matrix
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		ResidIter(R, Ac, bc, xc); //% residual correction
		C = A * R; //% right preconditioning
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				B(i, j) = (i == j) ? 1.0 - C(i, j) : -C(i, j);
			}
		}
		mag(reduce(B), Bm); //% |I-RA|
		if ((rho = rhoSpectral(Bm)) >= 1.0) { return NotStronglyRegular; }
		z = b - A * xc; //% residual correction
		y0 = z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (verification step)
		//%-------------------------------------------------------------
		while (true) {
			++niter;
			if (niter > MAXITER) { break; } //% maximal number of iterations exeeded
			blow(y0, eps, mi); //% blowing up is used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) {
					s = s + B(i, j) * y0[j];
				}
				y0[i] = z[i] + s;
			}
			if (reduce(y1).set_strictly_contains(reduce(y0))) {
				w = xc + R * y0;
				return Success; //% verified solution has been found
			}
		}
		if (dist1(y0, y1) < EPSX) {
			w = xc + R * y0;
			return NonVerifiedSolution;
		}
		w = xc + R * y0;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Rump <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using revised affine forms.
//% Linear dependencies are passed through iterations.
//% !!!Right!!! preconditioning
//%----------------------------------------------------------------------------
RumpRLCond(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	const dmatrix& L, // left preconditiong matrix
	const dmatrix& R, // right preconditioning matxix
	aafrvector& w, //% interval outer solution
	real &rho,
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows();
	aafrvector z{ n }, y{ n }, y1{ n };
	aafrmatrix B{ n, n }, C{ n, n };
	dvector bc{ n }, xc{ n };
	dmatrix Ac{ n, n }, Cm{ n, n };

	w = aafrvector{ n };
	mid(A, Ac); //% midpoint matrix
	mid(b, bc); //% midpoint vector
	gw_gausse(Ac, bc, xc); //% midpoint solution
	B = L * A * R; //% left and right preconditioning
	mid(B, Ac);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			C(i, j) = (i == j) ? 1.0 - B(i, j) : -B(i, j); //% C=I-RA
		}
	}
	mag(reduce(C), Cm); //% |I-RA|
	rho = rhoSpectral(Cm);
	if (rho >= 1.0) { return NotStronglyRegular; }
	z = L * (b - A * xc); //% residual correction and left-preconditioning
	y = z; //% starting point for the iteration
	niter = 0;
	//%-------------------------------------------------------------
	//% Main iteration loop (verification step)
	//%-------------------------------------------------------------
	while (true) {
		++niter;
		if (niter > MAXITER) { break; } //% maximal number of iterations exeeded
		blow(y, eps, mi); //% blowing up is used to obtain verified solution
		y1 = y;
		for (int i = 0; i < n; ++i) {
			AAFR s = 0.0;
			for (int j = 0; j < n; ++j) {
				s = s + C(i, j) * y[j];
			}
			y[i] = z[i] + s;
		}
		if (reduce(y1).set_strictly_contains(reduce(y))) {
			w = xc + R * y;
			return Success; //% verified solution has been found
		}
	}
	if (dist1(y, y1) < EPSX) {
		w = xc + R * y;
		return NonVerifiedSolution;
	}
	w = xc + R * y;
	return MaxIterExceeded;
} //% >>> Rump <<<

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% LU_decomposition (wihtout pivoting!!!!) to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
RumpRLCondLU(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int &niter)
{
	int n = A.num_rows(), res = Success;
	real d;
	cvector c{ n };
	dmatrix Ac{ n, n }, L{ n, n }, U{ n, n };
	mid(A, Ac);
	if (inv(Ac)) {
		LU_decompositionNP(Ac, c, d); //% LUDecomposition without pivoting
		LUComponents(Ac, L, U); //% extracting L and U from A
		res = RumpRLCond(A, b, U, L, w, rho, niter);
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% LU_decomposition with PIVOTING to obtain preconditioning matrices
//% Remark: Pivoting worsen the results
//%----------------------------------------------------------------------------
RumpRLCondLUP(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int &niter)
{
	int n = A.num_rows(), res = Success;
	real d;
	cvector c{ n };
	dmatrix Ac{ n, n }, L{ n, n }, U{ n, n }, P{ n, n }, Pinv{ n, n };
	aafrmatrix A2{ n, n };

	mid(A, Ac);
	if (inv(Ac)) {
		LU_decomposition(Ac, c, d); //% LUDecomposition with pivoting
		LUComponents(Ac, L, U); //% extracting L and U from A
		P.fill_in(0.0);
		for (int i = 0; i < n; ++i) {
			P(i, c[i]) = 1.0;
		}
		Pinv = P;
		inv(Pinv);
		A2 = A * Pinv;
		res = RumpRLCond(A2, b, U, L, w, rho, niter);
		w = Pinv * w;
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% SVD_decomposition to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
RumpRLCondSVD(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int & niter)
{
	int n = A.num_rows(), res = Success;
	cvector c{ n };
	dvector ww{ n };
	dmatrix Ac{ n, n }, L{ n, n }, U{ n, n }, V{ n, n };
	//
	mid(A, Ac);
	if (inv(Ac)) {
		SVDcmp(Ac, ww, V, n, n); //% singular value decomposition
		LUComponentsSVD(Ac, ww, V, L, U); //% extracting L and U from Ac, ww and V
		res = RumpRLCond(A, b, L, U, w, rho, niter);
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% QR_decomposition to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
RumpRLCondQR(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real & rho,
	int & niter)
{
	int n = A.num_cols(), res = Success;
	dvector c{ n }, d{ n };
	dmatrix R{ n, n }, L{ n, n }, U{ n, n };

	mid(A, R);
	if (inv(R)) {
		QRdcmp(R, n, n, c, d); //% QR decompotion
		LUComponentsQR(R, c, d, L, U); //% extracting L and U from R, c, and d
		res = RumpRLCond(A, b, U, L, w, rho, niter);
		return res;
	}
	return SingularMatrix;
}

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
	int& niter) //% number of iterations
{
	if (!oA[0].is_square()) return NonSquareMatrix;

	int	n = oA[0].num_rows();
	dvector	bc{ n }, xc{ n };
	dmatrix	Cm{ n, n }, R{ n, n };
	imatrix	A{ n, n }, B{ n, n }, C{ n, n };
	ivector	z{ n }, y0{ n }, y1{ n }, yk{ n };

	R = mid(oA, p); //% midpoint matrix
	bc = mid(ob, p); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		xc = R * bc; //% midpoint solution
		cpmatrix(R, p, oA, C); //% C = R*A(p)
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				B(i, j) = (i != j) ? -C(i, j) : 1.0 - C(i, j);
			}
		}
		mig(B, Cm); // Cm = <C>
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		zvector(xc, R, p, oA, ob, z); //% z = R*(b(p)-A(p)*x0)
		y0 = z; //% starting point for the iteration
		niter = 0;
		//% Main iteration loop
		while (true) {
			++niter;
			if (niter > MAXITER) { return Failed; }
			blow(y0, eps, mi);
			y1 = y0;
			for (int i = 0; i < n; ++i) {
				yk = z + (B * y1);
				y1[i] = yk[i];
			}
			if (y0.set_strictly_contains(yk)) {
				w = xc + yk;
				return Success;
			}
			y0 = yk;
		}
	}
	return SingularMatrix;
} //% >>> Rump with oA and ob <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump MRHS <<<
//% for systems with *** Multiple Right-hand Side (MHRS) ***
//% Self-verified fixed point iteration for solving parametric interval 
//% systems with multiple righ-hand side. The method is based on Krawczyk 
//% iteration.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& b, //% revised affine right-hand side matrix
	aafrmatrix& W, //% revised affine solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int	n = A.num_rows(), m = b.num_cols();
	dmatrix	Bc{ n, m }, Xc{ n, m };
	dmatrix	Cm{ n, n }, R{ n, n };
	aafrmatrix	B{ n, n }, C{ n, n };
	aafrmatrix	Z{ n, m }, Y0{ n, m }, Y1{ n, m }, yk{ n, m };
	aafrmatrix Ax{ n, m };
	AAFR ykk;

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% inverse of the midpoint matrix
		mid(b, Bc); //% midpoint vector
		Xc = R * Bc; //% midpoint solution
		B = R * A; //% left preconditioning
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				C(i, j) = (i != j) ? -B(i, j) : 1.0 - B(i, j);
			}
		}
		mag(reduce(C), Cm);
		if (rhoSpectral(Cm) >= 1.0) { return NotStronglyRegular; }
		Z = R * (b - A * Xc); //% right preconditioning and residual correction
		Y0 = Z; //% starting point for the iteration
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(Y0, eps, mi); //% blowing up is used to obtain verified solution
			Y1 = Y0;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					AAFR ykk = 0.0;
					for (int k = 0; k < n; ++k) {
						ykk = ykk + C(i, k) * Y0(k, j);
					}
					Y0(i, j) = Z(i, j) + ykk;
				}
			}
			if (StrictlyContains(reduce(Y1), reduce(Y0))) { //% solution found
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						W(i, j) = Xc(i, j) + Y1(i, j);
					}
				}
				return Success; //% verified solution has been found
			}
		} while (true);
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
	aafrvector w1{ n }, c{ n };
	aafrmatrix C{ n, n };
	dmatrix R{ n, n }, Br{ n, n }, R2{ n, n };

	mid(A, R); //% midpoint matrix
	//std::cout << "cond(A) = " << cond(R) << std::endl;
	if (inv(R)) {
		C = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		mid(C, R2);
		//std::cout << "condi(C) = " << cond(R2) << std::endl;
		c = R * b; //% O(n^2*K)
		rad(reduce(C), Br); //% R*A
		if (rhoSpectral(Br) >= 1.0) { 
			std::cout << rhoSpectral(Br) << std::endl;
			return NotStronglyRegular; 
		}
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < i; ++j) { //% lower triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; ++j) { //% upper triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
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
//% self verified method with pre-conditioning using the midpoint inverse
//% The iterative scheme of the method is as follows:
//% x=D^{-1}(b-(L+U)x)
//% It origines from the following decompsition of the system matrix:
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
//%---------------------------------------------------------------------------- 
GaussSeidel(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows();
	aafrvector w1{ n }, c{ n };
	aafrmatrix C{ n, n };
	dmatrix R{ n, n }, Br{ n, n }, Bc{ n, n }, Bcm{ n, n };

	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		c = b; //% O(n^2*K)
		rad(reduce(C), Br); //% |I-RA|
		mid(reduce(C), Bc);
		inv(Bc);
		mag(Bc, Bcm);
		Br = Bcm * Br;
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < i; ++j) { //% lower triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; ++j) { //% upper triangular
					if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + C(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
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
//% self verified method with pre-conditioning using the midpoint inverse
//% and residual correction
//% The iterative scheme of the method is as follows:
//% x=D^{-1}(b-(L+U)x), A = (L+D+U)
//% http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
//%---------------------------------------------------------------------------- 
ParamGaussSeidelRC(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }
	int n = A.num_rows();
	aafrvector w1{ n }, c{ n };
	aafrmatrix B{ n, n };
	dvector xc{ n }, bc{ n };
	dmatrix R{ n, n }, Br{ n, n }, Bc{ n, n }, Bcm{ n, n };

	mid(A, R); //% midpoint matrix
	mid(b, bc); //% midpoint of the right-hand vector
	if (inv(R)) { //% midpoint inverse
		xc = R * bc; //% midpoint solution
		B = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		rad(reduce(B), Br); //% rad(RA)
		mid(reduce(B), Bc);
		inv(Bc);
		mag(Bc, Bcm);
		Br = Bcm * Br;
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		c = R * (b - A * xc); //% O(n^2*K)
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < i; ++j) { //% lower triangular
					if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; ++j) { //% upper triangular
					if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				w[i] = (c[i] - s) / B(i, i); //% O(K*log(K))
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				w = xc + w;
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(w1, w) < EPSX) {
			w = xc + w;
			return NonVerifiedSolution;
		}
		w = xc + w;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> ParamGaussSeidel <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% self verified method with right pre-conditioning with midpoint inverse
//%----------------------------------------------------------------------------
ParamGaussSeidelRCond(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }

	int n = A.num_rows();
	aafrvector w1{ n };
	aafrmatrix B{ n, n };
	dmatrix	R{ n, n }, Br{ n, n };

	mid(A, R);
	if (inv(R)) {
		B = A * R; //% right pre-conditioning
		rad(reduce(B), Br); //% |I-RA|
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		w.fill_in(0.0);
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			for (int i = 0; i < n; ++i) { // all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
				AAFR s = 0.0;
				for (int j = 0; j < i; ++j) { //% lower triangular
					if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; ++j) { //% upper triangular
					if (B(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
						s = s + B(i, j) * w1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				w[i] = (b[i] - s) / B(i, i); //% O(K*log(K))
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
//% self verified method with pre-conditioning using the midpoint inverse
//% and residual correction
//% The iterative scheme of the method is as follows:
//% x=D^{-1}(b-(L+U)x), A = (L+D+U)
//% http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
//%---------------------------------------------------------------------------- 
ParamGaussSeidelRLCond(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	const dmatrix& L,
	const dmatrix& R,
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) { return NonSquareMatrix; }
	int n = A.num_rows();
	real rho;
	aafrvector w1{ n }, c{ n };
	aafrmatrix B{ n, n };
	dvector xc{ n }, bc{ n };
	dmatrix R0{ n, n }, Br{ n, n }, Bc{ n, n };

	mid(b, bc); //% midpoint of the right-hand vector
	mid(A, R0);
	if (inv(R0)) {
		xc = R0 * bc; //% midpoint solution
		B = L * A * R; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		mid(B, Bc);
		if (!inv(Bc)) { return SingularMatrix; }
		rad(reduce(B), Br); //% |I-RA|
		rho = rhoSpectral(Br);
		if (rho >= 1.0) { return NotStronglyRegular; }
		c = L * (b - A * xc); //% O(n^2*K)
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < i; ++j) { //% lower triangular
					if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; ++j) { //% upper triangular
					if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
						s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				w[i] = (c[i] - s) / B(i, i); //% O(K*log(K))
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				w = xc + R * w;
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(w1, w) < EPSX) {
			w = xc + R * w;
			return NonVerifiedSolution;
		}
		w = xc + R * w;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> ParamGaussSeidelLRCond <<<

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% LU_decomposition (wihtout pivoting!!!!) to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
PGSIRLCondLU(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int &niter)
{
	int n = A.num_rows(), res = Success;
	real d;
	cvector c{ n };
	dmatrix Ac{ n, n }, L{ n, n }, U{ n, n };
	mid(A, Ac);
	if (inv(Ac)) {
		LU_decompositionNP(Ac, c, d); //% LUDecomposition without pivoting
		LUComponents(Ac, L, U); //% extracting L and U from A
		res = ParamGaussSeidelRLCond(A, b, U, L, w, niter);
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% LU_decomposition (wihtout pivoting!!!!) to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
PGSIRLCondSVD(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int &niter)
{
	int n = A.num_rows(), res = Success;
	cvector c{ n };
	dvector ww{ n };
	dmatrix Ac{ n, n }, L{ n, n }, U{ n, n }, Acc{ n,n }, V{ n, n };
	mid(A, Ac);
	Acc = Ac;
	if (inv(Ac)) {
		SVDcmp(Ac, ww, V, n, n); //% singular value decomposition
		LUComponentsSVD(Ac, ww, V, L, U); //% extracting L and U from Ac, ww and V
		res = ParamGaussSeidelRLCond(A, b, L, U, w, niter);
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Krawczyk iteration with double preconditioning
//% The inverse of the midpoint matrix is decomposed by using
//% LU_decomposition (wihtout pivoting!!!!) to obtain preconditioning matrices
//% Pivoting worsen the resuls
//%----------------------------------------------------------------------------
PGSIRLCondQR(const aafrmatrix& A,
	const aafrvector& b,
	aafrvector& w,
	real &rho,
	int &niter)
{
	int n = A.num_rows(), res = Success;
	dvector ww{ n }, c{ n }, d{ n };
	dmatrix R{ n, n }, L{ n, n }, U{ n, n }, Acc{ n,n }, V{ n, n };
	mid(A, R);
	if (inv(R)) {
		QRdcmp(R, n, n, c, d); //% QR decompotion
		LUComponentsQR(R, c, d, L, U); //% extracting L and U from R, c, and d
		res = ParamGaussSeidelRLCond(A, b, U, L, w, niter);
		return res;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel MRHS <<<
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
	aafrmatrix w1{ n, m }, c{ n, m };
	aafrmatrix B{ n, n };
	dmatrix R{ n, n }, Br{ n, n };

	w = aafrmatrix{ n, m };
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		B = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
		rad(reduce(B), Br); //% |I-RA|
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		c = R * b; //% O(n^2*K)
		w = c;
		niter = 0;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(w, eps, mi); //% blowing up is used to obtain verified solution
			w1 = w;
			//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < m; ++k) {
					AAFR s = 0.0;
					for (int j = 0; j < i; ++j) { //% lower triangular
						if (B(i, j).reduce() != 0.0 && w(j, k).reduce() != 0.0) {
							s = s + B(i, j) * w(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					for (int j = i + 1; j < n; ++j) { //% upper triangular
						if (B(i, j).reduce() != 0.0 && w1(j, k).reduce() != 0.0) {
							s = s + B(i, j) * w(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					w(i, k) = (c(i, k) - s) / B(i, i); //% O(K*log(K))
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
	aafrvector c{ n }, v{ n }, v1{ n };
	aafrmatrix B{ n, n }, Ap{ n, n };
	dmatrix R{ n, n }, Br{ n, n }, Ac{ n, n };

	mid(A, Ac);
	R = Ac;
	if (inv(R)) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Ap(i, j) = A(i, j) - Ac(i, j);
			}
		}
		//% left pre-conditioning
		B = R * Ap;
		rad(reduce(B), Br); //% |I-RA|
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		c = R * b;
		v = c;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) { break; }
			blow(v, eps, mi); //% blowing up is used to obtain verified solution
			v1 = v;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) {
					if (B(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + B(i, j) * v[j];
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
	aafrvector d{ n }, v{ n }, v1{ n };
	aafrmatrix Ap{ n, n }, B{ n, n };
	dvector bc{ n }, xc{ n };
	dmatrix Br{ n, n }, R{ n, n }, Ac{ n, n };

	mid(b, bc); //% midpoint vector mid(b)
	mid(A, Ac);
	R = Ac;
	if (inv(R)) {
		xc = R * bc; //% midpoint solution
		d = R * (b - A * xc); //% residual correction
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Ap(i, j) = A(i, j) - Ac(i, j);
			}
		}
		B = R * Ap;
		mag(reduce(B), Br); //% |I-RA|
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		//% transform initial interval enclosure into revised affine form
		v = d;
		//%-------------------------------------------------------------
		//% Main iteration loop (+verification step)
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) break;
			blow(v, eps, mi);
			v1 = v;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) { //% use already updated entries
					if (B(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + B(i, j) * v[j];
					}
				}
				v[i] = d[i] - s;
			}
			if (reduce(v1).set_strictly_contains(reduce(v))) {
				w = xc + v;
				return Success; //% verified solution has been found
			}

		} while (true);
		if (dist1(v1, v) < EPSX) {
			w = xc + v;
			return NonVerifiedSolution;
		}
		w = xc + v;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> MKolevRC <<<

int
//%----------------------------------------------------------------------------
//% Modified version of Kolev's method
//% computing with reduced affine forms
//% with !!!right!!! pre-conditioning
//%----------------------------------------------------------------------------
MKolevRCond(const aafrmatrix& A, //% revied affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w, //% p-solution (revised affine matrix)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector d(n), v(n), v1(n);
	aafrmatrix Ap(n, n), B(n, n);
	cvector	xpidx(1);
	dvector x0(n), xpcfs(1);
	dmatrix R(n, n), Am{ n, n }, Br{ n, n };

	mid(A, Am);
	R = Am;
	if (inv(R)) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Ap(i, j) = A(i, j) - Am(i, j);
			}
		}
		B = Ap * R; //% right pre-conditioning of system matrix
		rad(reduce(B), Br); //% |I-RA|
		if (rhoSpectral(Br) >= 1.0) { return NotStronglyRegular; }
		v = b;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) break;
			blow(v, eps, mi);
			v1 = v;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) {
					if (B(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + B(i, j) * v[j];
					}
				}
				v[i] = b[i] - s;
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
MKolevRCondRC(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations 
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector d{ n }, v{ n }, v1{ n };
	dvector bm{ n }, xc{ n };
	aafrmatrix B{ n, n }, Ap{ n, n };
	dmatrix	Am{ n, n }, R{ n, n };

	mid(A, Am);
	R = Am;
	if (inv(R)) {
		mid(b, bm); //% midpoint vector
		xc = R * bm; //% midpoint solution
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Ap(i, j) = A(i, j) - Am(i, j);
			}
		}
		B = Ap * R; //% right pre-conditioning
		d = (b - A * xc); //% residual correction
		v = d;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) break;
			blow(v, eps, mi);
			v1 = v;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) { //% use already updated entries
					if (B(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
						s = s + B(i, j) * v[j];
					}
				}
				v[i] = d[i] - s;
			}
			if (reduce(v1).set_strictly_contains(reduce(v))) {
				w = xc + R * v;
				return Success;
			}
		} while (true);
		if (dist1(v1, v) <= EPSX) {
			w = xc + R * v;
			return NonVerifiedSolution;
		}
		w = xc + R * v;
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> MKolevRCII <<<

int
//%----------------------------------------------------------------------------
//% >>> Jacobi iteration <<<
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%---------------------------------------------------------------------------- 
Jacobi(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& w,	//% p-solution (revised affine vector)
	int& niter)	//% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector w1{ n }, c{ n };
	aafrmatrix C{ n, n };
	dmatrix R{ n, n };

	mid(A, R);
	if (inv(R)) {
		C = R * A;
		c = R * b;
		w = c;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) break;
			blow(w, eps, mi);
			w1 = w;
			for (int i = 0; i < n; ++i) {
				AAFR s = 0.0;
				for (int j = 0; j < n; ++j) {
					if (i != j) {
						if (C(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
							s = s + C(i, j) * w[j];
						}
					}
				}
				w[i] = (c[i] - s) / C(i, i);
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				return Success; //% verified solution has been found
			}
		} while (true);
		if (dist1(w, w1) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> Jacobi <<<

int
//%----------------------------------------------------------------------------
//% Succesive overrelaxation method
//% Iterative method (preconditioning with minpoint inverse)
//%---------------------------------------------------------------------------- 
SOR(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand vector
	aafrvector& w, //% interval solution
	real wgt, //% relaxation parameter
	int& niter)	//% number of iterations
{
	int n = A.num_rows();
	ivector xp(n);
	dmatrix R(n, n);
	aafrvector w1(n), c(n);
	aafrmatrix C(n, n), Rf(n, n);

	mid(A, R);
	if (inv(R)) {
		C = R * A;
		c = R * b;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			++niter;
			if (niter > MAXITER) break;
			blow(w, eps, mi);
			w1 = w;
			for (int i = 0; i < n; ++i) {
				AAFR s1 = 0.0;
				AAFR s2 = 0.0;
				for (int j = 0; j < i; ++j) {
					s1 = s1 + C(i, j) * w[j];
				}
				for (int j = i + 1; j < n; ++j) {
					s2 = s2 + C(i, j) * w[j];
				}
				w[i] = (1.0 - wgt)*w[i] + wgt * (c[i] - s1 - s2) / C(i, i);
			}
			if (reduce(w1).set_strictly_contains(reduce(w))) {
				return Success;
			}
		} while (true);
		if (dist1(w, w1) < EPSX) {
			return NonVerifiedSolution;
		}
		return MaxIterExceeded;
	}
	return SingularMatrix;
} //% >>> SuccessiveOverrelaxation <<<
