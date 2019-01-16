//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Self-verified fixed point iteration
//%
//%     Contains: 
//%     CheckForZeros, Accurate, StrictlyContains, ResidIter
//%     Rump for PILS, Rump for MRHS, Rump for ILS
//%
//%##########################################################################
//%

#include <iostream>
#include <math.h>
#include <time.h>
#include <conio.h>
#include "../utils/stdafx.h"
#include "../interval/interval_base.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/randnum.h"
#include "../utils/inverse.h"
#include "../truss/TrussStructure.h"
#include "../affine/aafr.h"
#include "rump.h"

const real PGSIeps = 1.0e-8;
const real eps = 1.0e-9;
const real mi = spn(); //% smallest positive number

#define SIGN(x) (((x) > 0) ? 1 : ((x) < 0) ? -1 : 0)
#define MAXITERMP 100

static void
//%----------------------------------------------------------------------------
//% The vectors x and y are successive approximations for the solution of a
//% linear system of equations computed by iterative refinement. If a compo-
//% nent of y is diminished by more than 'Factor', it is a good candidate for
//% a zero entry. Thus, it is set to zero.
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

bool
//%----------------------------------------------------------------------------
//% >>> StrictlyContains <<<
//% Verifies if interval matrix [A] is strictly contained in matrix [B] 
//%----------------------------------------------------------------------------
StrictlyContains(const imatrix& A, const imatrix& B)
{
	assert(B.num_cols() == A.num_cols() && B.num_rows() == A.num_rows());

	int n = A.num_rows(), m = A.num_cols();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (!B(i, j).set_strictly_contains(A(i, j))) {
				return false;
			}
		}
	}
	return true;
} //% >>> StrictlyContains <<<

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
		if (ok || k > MAXITERMP) return;
		x = y;
		k++;
	} while (true);
} //% >>> ResidIter <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using revised affine forms.
//% Linear dependencies are passed through iterations.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% interval outer solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	aafrvector z(n), y0(n), y1(n);
	aafrmatrix B(n, n), C(n, n); //% this is the main difference
	dvector bc(n), x0(n);
	dmatrix Ac(n, n), Cm(n, n), R(n, n), I(n, n);

	w = aafrvector(n);
	mid(A, Ac); //% midpoint matrix
	R = Ac;
	if (inv(R)) {
		mid(b, bc); //% midpoint vector
		x0 = R * bc; //% midpoint solution
		ResidIter(R, Ac, bc, x0);
		C = R * A; //% preconditioning
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - C(i, j);
			}
		}
		mig(reduce(B), Cm);
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		z = R * (b - A * x0);
		y0 = z; //% starting point for the iteration
		niter = 0;
		//% Main iteration loop (verification step)
		while (true) {
			niter++;
			if (niter > MAXITERMP) return Failed;
			blow(y0, eps, mi); //% used to obtain verified solution
			y1 = y0;
			for (int i = 0; i < n; i++) {
				AAFR ykk = 0.0;
				for (int j = 0; j < n; j++) {
					ykk = ykk + B(i, j) * y1[j];
				}
				y1[i] = z[i] + ykk;
			}
			if (reduce(y0).set_strictly_contains(reduce(y1))) { //% solution found
				for (int i = 0; i < n; i++) {
					w[i] = x0[i] + y1[i];
				}
				return Success;
			}
			y0 = y1;
		}
	}
	return SingularMatrix;
} //% >>> Rump <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump for systems with Multiple Right-hand Side (MHRS) <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems with multiple righ-hand side. The method is based on Krawczyk 
//% iteration.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& b, //% revised affine right-hand side matrix
	aafrmatrix& w, //% revised affine solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), m = b.num_cols();
	dmatrix	Bc(n, m), X0(n, m);
	dmatrix	Cm(n, n), R(n, n), I(n, n);
	aafrmatrix	B(n, n), C(n, n);
	aafrmatrix	Z(n, m), Y0(n, m), Y1(n, m), yk(n, m), yy(n, m), D(n, m);
	aafrmatrix Ax(n, m);
	AAFR ykk;

	mid(A, R);
	if (inv(R)) {
		mid(b, Bc); //% midpoint vector
		X0 = R * Bc; //% midpoint solution
		C = R * A; //% preconditioning
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - C(i, j);
			}
		}
		mig(reduce(B), Cm);
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		Z = R * (b - A * X0);
		Y0 = Z; //% starting point for the iteration
		niter = 0;
		//% Main iteration loop (verification step)
		while (true) {
			niter++;
			if (niter > MAXITERMP) return Failed;
			blow(Y0, eps, mi); //% used to obtain verified solution
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
				return Success;
			}
			Y0 = Y1;
		}
	}
	return SingularMatrix;
} //% >>> Rump for MHRS <<<

int
//%---------------------------------------------------------------------------- 
//% >>> Rump for ILS <<<
//% Verified (Rump) fixed-point iteration for solving interval linear
//% systems (ILS).
//%----------------------------------------------------------------------------
Rump(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand side vector
	ivector& x) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), stepCount = 0;
	dvector	bc(n), x0(n);
	dmatrix	R(n, n), I(n, n), Bm(n, n);
	ivector	z(n), y0(n), yk(n), yy(n);
	imatrix	B(n, n), C(n, n);

	mid(A, R); //% midpoint matrix
	x = ivector(n);
	if (inv(R)) { //% midpoint inverse
		mid(b, bc);	//% midpoint vector
		x0 = R * bc; //% midpoint solution
		Unit(I);
		B = I - R * A;
		mig(B, Bm); //% Ostrowski comparison matrix
		if (rhoSpectral(Bm) >= 1.0) return NotStronglyRegular;
		z = R * (b - A * x0); //% residual correction
		y0 = z;
		while (true) {
			stepCount++;
			blow(y0, eps, mi);
			yk = z + (B * y0);
			if (y0.set_strictly_contains(yk)) {
				x = x0 + yk;
				return Success;
			}
			y0 = yk;
		}
	}
	return SingularMatrix;
} //% >>> Rump for ILS <<<
