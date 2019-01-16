//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: 
//%     Iterative methods for solving interval interval linear systems
//%
//%     Contains:
//%     firstSOlution
//%     Rump for ILS
//%     GaussSeidel, GaussSeidel RC
//%     JacobiCond, Jacobi, JacobiRC,
//%     Krawczyk, KrawczykRC,
//%
//%##########################################################################
//%

#include <conio.h>
#include "../utils/stdafx.h"
#include "../interval/interval_base.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../param_solver/direct_solvers.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/iterative_isolvers.h"


const real PGSIeps = 1.0e-8;
const real eps = 1.0e-9;
const real mi = spn(); //% smallest positive number

#define SIGN(x)		(((x) > 0) ? 1 : ((x) < 0) ? -1 : 0)
#define EPSX		1.0e-11 //% accuracy
#define MAXITER		200 //% maximal number of iterations

int
//%----------------------------------------------------------------------------
//% >>> firstSOlution <<<
//% Initital solution used by iterative methods for solving 
//% interval linear systems. Method described in Neumaier's book.
//%----------------------------------------------------------------------------
firstSolution(const dmatrix& R, //% real matrix (preconditioner)
	const imatrix& A, //% interval system matrix
	const ivector& b, //% interval vector
	ivector& w)	//% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size();
	real alpha;
	dvector u(n), eps(n), bmag(n);
	dmatrix Cm(n, n), B(n, n);

	w.fill_in(0.0);
	eps.fill_in(0.1e-15);
	mig(A, Cm);
	B = Cm;
	if (inv(Cm)) {
		mag(b, bmag);
		u = Cm * (bmag + eps);
		alpha = 0.01;
		if (B*u - alpha * bmag < 0.0) {
			return SingularMatrix;
		}
		w = (1.0 / alpha) * u * interval(-1.0, 1.0);
		return Success;
	}
	return SingularMatrix;
} //% >>> firstSOlution <<<

int
//%----------------------------------------------------------------------------
//% >>> Krawczyk iteration <<<
//% for solving interval linear systems
//%----------------------------------------------------------------------------
Krawczyk(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w, //% interval solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	real eps = 1.0e-5, mi = spn(), rs;
	dvector bm(n), x0(n);
	dmatrix	R(n, n), I(n, n), Br(n, n);
	ivector	w0(n), c(n), y0(n), y1(n);
	imatrix C(n, n), B(n, n);

	w = ivector(n);
	mid(b, bm);	//% midpoint vector
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		x0 = R * bm; //% midpoint solution
		C = R * A; //% left pre-conditioning
		c = R * b; //% left pre-conditioning
		Unit(I);
		B = I - C;
		mag(B, Br);
		rs = rhoSpectral(Br);
		if (rs >= 1.0) return NotStronglyRegular;
		y0 = c;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			y1 = y0;
			for (int i = 0; i < n; i++) {
				interval s = 0.0;
				for (int j = 0; j < n; j++) {
					s = s + B(i, j) * y0[j];
				}
				y0[i] = c[i] + s;
			}
		} while (distH(y0, y1) > EPSX);
		w = y0;
		return Success;
	}
	return SingularMatrix;
} //% >>> Krawczyk <<<

int
//%----------------------------------------------------------------------------
//% >>> Krawczyk iteration <<<
//% with residual correction
//% for solving interval linear systems
//%----------------------------------------------------------------------------
KrawczykRC(const imatrix&A, //% revised affine system matrix
	const ivector& b, //% revised affine right-hand vector
	ivector& w, //% interval solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector bm(n), x0(n);
	dmatrix	R(n, n), I(n, n);
	ivector	w0(n);
	imatrix C(n, n), B(n, n);
	ivector c(n), y0(n), y1(n);

	mid(b, bm);	//% midpoint vector
	mid(A, R); //% midpoint matrix
	w = ivector(n);
	if (inv(R)) { //% midpoint inverse
		x0 = R * bm; //% midpoint solution
		C = R * A; //% left pre-conditioning
		c = R * (b - A * x0); //% left pre-conditioning and residual correction
		Unit(I);
		B = I - C;
		y0 = c; // z = R * (b - A * x0); //% residual correction
				//%-------------------------------------------------------------
				//% Main iteration loop
				//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			y1 = y0;
			for (int i = 0; i < n; i++) {
				interval s = 0.0;
				for (int j = 0; j < n; j++) {
					s = s + B(i, j) * y0[j];
				}
				y0[i] = c[i] + s;
			}
		} while (distH(y0, y1) > EPSX);
		w = x0 + y0;
		return Success;
	}
	return SingularMatrix;
} //% >>> KrawczykRC <<<

real
//%----------------------------------------------------------------------------
//% >>> Jacobi condition <<<
//% condition for convergence of Jacobi method
//%----------------------------------------------------------------------------
JacobiCond(const imatrix& C)
{
	int n = C.num_cols();
	dmatrix Cr(n, n), Am(n, n);
	imatrix Dd(n, n), Ap(n, n);
	Dd.fill_in(0.0);
	Ap.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) Dd(i, j) = C(i, j); else Ap(i, j) = C(i, j);
		}
	}
	mag(Dd, Am);
	Cr = Am;
	mag(Ap, Am);
	Cr = Cr * Am;
	return rhoSpectral(Cr);
}

int
//%----------------------------------------------------------------------------
//% >>> Jacobi iteration <<<
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%----------------------------------------------------------------------------
Jacobi(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& x0, //% interval solution vector
	int& niter) //% number of iterations 
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	real rs;
	ivector	x1(n), D(n), bb(n);
	dvector	xinf(n), xsup(n), xmid(n);
	ivector	xp(n), wek(n);
	imatrix	C(n, n);
	dmatrix	R(n, n);

	//% Initial enclosure computed using ISM
	mid(A, R);
	if (inv(R)) { //% midpoint inverse
		C = R * A; //% left pre-conditioning
		rs = JacobiCond(C);
		if (rs >= 1.0) return NotStronglyRegular;
		bb = R * b; //% left pre-conditioning
		err = firstSolution(R, C, bb, wek); //% initial solution
		if (err != Success) return err;
		x0 = wek;
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			x1 = x0;
			for (int i = 0; i < n; i++) {
				interval s = 0.0;
				for (int j = 0; j < n; j++) {
					if (i != j) {
						if (C(i, j) != 0.0 && x1[j] != 0.0) {
							s = s + C(i, j) * x1[j];
						}
					}
				}
				x0[i] = (bb[i] - s) / C(i, i);
			}
		} while (dist1(x0, x1) > EPSX);
		return Success;
	}
	return SingularMatrix;
} //% >>> Jacobi iteration <<<

int
//%----------------------------------------------------------------------------
//% >>> Jacobi iteration <<<
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%---------------------------------------------------------------------------- 
Jacobi(const dmatrix& R, //% conditioning matrix
	const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& x0, //% interval solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	real rs;
	ivector	x1(n), D(n), bb(n);
	dvector	xinf(n), xsup(n);
	ivector	xp(n), wek(n);
	imatrix	C(n, n);

	//% Initial enclosure computed using ISM
	C = R * A; //% left pre-conditioning
	rs = JacobiCond(C);
	if (rs >= 1.0) return NotStronglyRegular;
	bb = R * b; //% left pre-conditioning
	err = firstSolution(R, C, bb, wek); //% initial solution
	if (err != Success) return err;
	x0 = wek;
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		for (int i = 0; i < n; i++) {
			interval s = 0.0;
			for (int j = 0; j < n; j++) {
				if (i != j) {
					if (C(i, j) != 0.0 && x1[j] != 0.0) {
						s = s + C(i, j) * x1[j];
					}
				}
			}
			x0[i] = (bb[i] - s) / C(i, i);
		}
	} while (dist1(x0, x1) > EPSX);
	return Success;
} //% >>> Jacobi iteration <<<

int
//%----------------------------------------------------------------------------
//% Jacobi iteration with residual correction
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%----------------------------------------------------------------------------
JacobiRC(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& x0, //% interval solution vector
	int& niter) //% number of iterations 
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	ivector bb(n);
	dvector	bm(n), xc(n);
	dmatrix	R(n, n);

	mid(A, R);
	mid(b, bm);
	x0 = ivector(n);
	if (inv(R)) {
		xc = R * bm;
	//	Noise(xc);
		bb = b - A * xc; //% residual correction
		err = Jacobi(R, A, bb, x0, niter);
		if (err != Success) return err;
		x0 = xc + x0;
		return Success; //% success!!!
	}
	return SingularMatrix; //% singular matrix
} //% >>> Jacobi resid <<<

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

int
//%----------------------------------------------------------------------------
//% Jacobi iteration for standard interval linear systems
//% without taking into account dependencies 
//% TRZEBA PRZEROBIC NA WERSJE "ELEMENTOWA"
//%----------------------------------------------------------------------------
Jacobi(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	real error;
	int n = A.num_rows(), err;
	dvector xinf(n), xsup(n);
	dmatrix R(n, n);
	ivector x0(n), x1(n);
	ivector D(n), bb(n);
	imatrix L(n, n), U(n, n), Dinv(n, n), AA(n, n);

	err = firstSolution(R, AA, bb, x0); //% initial solution
	if (err != Success) return err;
	//x0.fill_in(0.0); //% initial point
	L.fill_in(0.0); //% lower triangular
	U.fill_in(0.0); //% upper triangular
	Dinv.fill_in(0.0);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		AA = R * A; //% left pre-conditioning
		bb = R * b; //% left pre-conditioning
		for (int i = 0; i < n; i++) {
			D[i] = AA(i, i);
			Dinv(i, i) = 1.0 / D[i];
			for (int j = 0; j < i; j++) {
				L(i, j) = AA(i, j);
				U(j, i) = AA(j, i);
			}
		}
		//%-------------------------------------------------------------
		//% Main iteration loop
		//%-------------------------------------------------------------
		do {
			error = 0.0;
			x1 = x0;
			x0 = Dinv * (bb - (L + U)*x0);
			for (int i = 0; i < n; i++) {
				error += ISLIB_ABS(x0[i].inf() - x1[i].inf()) + ISLIB_ABS(x0[i].sup() - x1[i].sup());
			}
		} while (error > EPSX);
		w = x0;
		return Success;
	}
	return SingularMatrix;
} //% >>> Jacobi <<<

int
//%----------------------------------------------------------------------------
//% >>> Preconditioned Gauss-Seidel iteration <<<
//% without residual correction
//%----------------------------------------------------------------------------
GaussSeidel(const imatrix& A, //% interval matrix
	const ivector& b, //% interval vector
	ivector& wek, //% interval solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), err;
	interval tmpwek, ttmp;
	dmatrix R(n, n);
	ivector tmp(n), bb(n);
	imatrix AA(n, n);

	mid(A, R); //% midpoint matrix = preconditioner
	inv(R); //% midpoint inverse
	AA = R * A; //% left pre-conditioning
	bb = R * b; //% left pre-conditioning
	err = firstSolution(R, AA, bb, wek); //% initial solution
	if (err != Success) return err;
	tmp = wek;
	niter = 0;
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	do {
		niter++;
		for (int i = 0; i < n; i++) {
			interval s = 0.0;
			for (int j = 0; j < n; j++) {
				if (i != j) {
					s += AA(i, j) * wek[j];
				}
			}
			ttmp = (bb[i] - s) / AA(i, i);
			ttmp.intersects(tmp[i], tmpwek);
			if (tmpwek.is_empty()) {
				return Failed;
			}
			wek[i] = tmpwek;
			real ss = 0.0;
			for (int k = 0; k < n; k++) {
				ss += ISLIB_ABS(tmp[k].inf() - wek[k].inf()) + ISLIB_ABS(tmp[k].sup() - wek[k].sup());
			}
			if (ss < EPSX) { //% stopping criterion
				return Success;
			}
			tmp[i] = wek[i];
		}
	} while (true);
} //% >>> GaussSeidel <<<

int
//%----------------------------------------------------------------------------
//% Pre-conditioned Gauss-Seidel iteration with residual correction
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%---------------------------------------------------------------------------- 
GaussSeidelRC(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& x0, //% interval solution vector
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	ivector bb(n);
	dvector	bm(n), xc(n);
	dmatrix	R(n, n);

	mid(A, R);
	mid(b, bm);
	x0 = ivector(n);
	if (inv(R)) {
		xc = R * bm;
		bb = b - A * xc; //% residual correction
		err = GaussSeidel(A, bb, x0, niter);
		if (err != Success) return err;
		x0 = xc + x0;
		return Success;
	}
	return SingularMatrix;
} //% >>> GaussSeidel resid <<<