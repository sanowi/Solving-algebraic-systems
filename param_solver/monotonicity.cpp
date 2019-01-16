//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Monotonicity Approach methods based on testing the sign
//%             of partial derivatives on an interval.
//%
//%##########################################################################
//%

#include <iostream>
#include <conio.h>
#include <omp.h>
#include <ctime>
#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/gausselim.h"
#include "../utils/inverse.h"
#include "../newton/Newton.h"
#include "../affine/aaf.h"
#include "../affine/aafr.h"
#include "../param_solver/kolev.h"
#include "../param_solver/degrauwe.h"
#include "../param_solver/iterative_psolvers.h"
#include "../param_solver/rump.h"
#include "../param_solver/old_methods.h"
#include "../examples/pexamples.h"

void
//%----------------------------------------------------------------------------
//% >>> PartAAFR <<<
//% Computes partial derivatives d_A(p)/d_pk in a system matrix
//% there are only linear dependencies, so the derivatives are constant values
//%----------------------------------------------------------------------------
PartAAFR(const aafrmatrix& Ar, //% revised affine system matrix
	dmatrix& A, //% real matrix (of derivatives)
	int dk) //% number of variable 
{
	int n = Ar.num_rows();
	cvector idx;
	dvector cfs;

	A.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			idx = Ar(i, j).idx();
			cfs = Ar(i, j).cfs();
			int K = idx.size();
			if (dk <= idx[K - 1]) {
				for (int kk = 1; kk < K; kk++) {
					if (idx[kk] == dk) {
						A(i, j) = cfs[kk];
					}
				}
			}
		}
	}
} //% >>> PartAAFR <<<

void
//%----------------------------------------------------------------------------
//% >>> PartAAFR <<<
//% Computes partial derivatives d_b(p)/d_pk in a right-hand vector
//% there are only linear dependencies, so the derivatives are constant values
//%----------------------------------------------------------------------------
PartAAFR(const aafrvector& br, //% revised affine vector
	dvector& b, //% real vector (of derivatives)
	int dk) //% number of variable 
{
	int n = br.size();
	cvector idx;
	dvector cfs;

	b.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		idx = br[i].idx();
		cfs = br[i].cfs();
		int K = idx.size();
		if (dk <= idx[K - 1]) {
			for (int kk = 1; kk < K; kk++) {
				if (idx[kk] == dk) {
					b[i] = cfs[kk];
				}
			}
		}
	}
} //% >>> PartAAFR <<<

void
//%----------------------------------------------------------------------------
//% Reduces affine matrix after computing derivatives
//%----------------------------------------------------------------------------
MonotReduce(int m, //% 0 - minimum, 1 - maximum
	const aafrmatrix& Ar, //% revised affine matrix
	const dvector dx, //% real vector
	aafrmatrix& A) //% revised affine matrix after reduction
{
	int n = Ar.num_cols();
	int K = dx.size();
	dvector cfs;
	cvector idx;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			idx = Ar(i, j).idx();
			cfs = Ar(i, j).cfs();
			int KK = idx.size();
			//% loop over the elements of affine form Ar(i, j)
			for (int k = 1; k < KK; k++) {
				if (dx[idx[k] - 1] == -1.0) {
					cfs[0] += (m == 1) ? -cfs[k] : cfs[k];
					cfs[k] = 0.0;
				}
				if (dx[idx[k] - 1] == 1.0) {
					cfs[0] += (m == 1) ? cfs[k] : -cfs[k];
					cfs[k] = 0.0;
				}
			}
			A(i, j) = AAFR(idx, cfs, Ar(i, j).r());
			A(i, j).reduceZeroes();
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Reduces affine vector after computing derivatives
//%----------------------------------------------------------------------------
MonotReduce(int m, //% 0 - minimu, 1 - maximum
	const aafrvector& br, //% revised affine vector
	const dvector dx, //% real vector
	aafrvector& b)  //% revised affine vector after reduction
{
	int n = br.size();
	int K = dx.size();
	dvector cfs;
	cvector idx;

	for (int i = 0; i < n; i++) {
		idx = br[i].idx();
		cfs = br[i].cfs();
		int KK = idx.size();
		//% loop over the elements of affine form br(i, j)
		for (int k = 1; k < KK; k++) {
			if (dx[idx[k] - 1] == -1) {
				cfs[0] += (m == 1) ? -cfs[k] : cfs[k];
				cfs[k] = 0.0;
			}
			if (dx[idx[k] - 1] == 1) {
				cfs[0] += (m == 1) ? cfs[k] : -cfs[k];
				cfs[k] = 0.0;
			}
		}
		b[i] = AAFR(idx, cfs, br[i].r());
		b[i].reduceZeroes();
	}
}

int
//%----------------------------------------------------------------------------
//% >>> Standard Monotonicity Approach <<<
//% The first argument of the method is a pointer to an iterative method that
//% is used to compute the derivatives on an interval.
//%----------------------------------------------------------------------------
MonotStandardI(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), //% function used to solve PILS
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& ww) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), precond = 0, err;
	int K = AAFR::getlast(); //% number of initial parameters e_1, ..., e_K
	int niter = 0;
	aafrvector xa(n), c(n), x0(n), xf(n), xff(n), bi(n);
	dvector pb(n);
	ivector w(n), wb(n);
	ovector dx(n);
	aafrmatrix Ai(n, n);
	dmatrix pA(n, n);

	ww = ivector(n);
	err = f(A, b, xa, niter); //% xa = x(p) p-solution to the system
	if (err != Success && err != NonVerifiedSolution) return err;
	//% compute derivatives over p_k
	wb = reduce(xa); //% xb = [x(p)] - interval solution of the system
	//% transform interval solution into revised affine vector
	xf.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		xf[i] = xf[i] + wb[i];
	}
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	//% compute derivatives over p_k
	//% compute right-hand vector of system with derivatives
	for (int k = 1; k < K + 1; k++) {
		PartAAFR(A, pA, k);
		PartAAFR(b, pb, k);
		c = pb - pA * xf; // xf; - here is the problem
		err = f(A, c, x0, niter);
		if (err != Success && err != NonVerifiedSolution) return err;
		w = reduce(x0);
		for (int i = 0; i < n; i++) {
			if (w[i] > 0.0) {
				dx[i][k - 1] = 1.0;
			}
			else if (w[i] < 0.0) {
				dx[i][k - 1] = -1.0;
			}
			else {
				dx[i][k - 1] = 0.0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		//% minimum
		MonotReduce(-1, A, dx[i], Ai);
		MonotReduce(-1, b, dx[i], bi);
		err = f(Ai, bi, xf, niter);
		if (err != Success && err != NonVerifiedSolution) return err;
		//% maximum
		MonotReduce(1, A, dx[i], Ai);
		MonotReduce(1, b, dx[i], bi);
		err = f(Ai, bi, xff, niter);
		if(err != Success && err != NonVerifiedSolution) return err;
		ww[i] = interval(xf[i].reduce().inf(), xff[i].reduce().sup());
	}
	return Success;
} //% >>> Standard Monotonicity Approach <<<

int
//%----------------------------------------------------------------------------
//% >>> Monotonicity approach <<<
//% based on augmented system matrix.
//% Gives the best results among all MAs, but is expensive.
//%----------------------------------------------------------------------------
MonotAugmentedD(int(*f)(const aafrmatrix&, const aafrvector&, ivector&), //% function used to solve PILS
	const aafrmatrix& A,  //% revised affine system matrix
	const aafrvector& b,  //%revised affine right-hand vector
	ivector& ww) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), precond = 0, err;
	int K = AAFR::getlast(); //% number of initial parameters e_1, ..., e_K
	int niter = 0;
	aafrvector xa(n), bi(n), bg(2 * n);
	aafrmatrix Ai(n, n), Ag(2 * n, 2 * n);
	dvector pb(n);
	ivector w(2 * n), xi(2 * n), xii(2 * n);
	ovector dx(n);
	dmatrix pA(n, n);

	ww = ivector(n);
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	Ag.fill_in(0.0);
	bg.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		bg[i] = b[i];
		for (int j = 0; j < n; j++) {
			Ag(i, j) = A(i, j);
			Ag(i + n, j + n) = A(i, j);
		}
	}
	//% computing derivatives over p_k
	for (int k = 1; k < K + 1; k++) {
		//% computing right-hand vector of system with derivatives
		PartAAFR(A, pA, k);
		PartAAFR(b, pb, k);
		for (int i = 0; i < n; i++) {
			bg[i + n] = pb[i];
			for (int j = 0; j < n; j++) {
				Ag(i + n, j) = pA(i, j);
			}
		}
		err = f(Ag, bg, w);
		if (err != Success) return err;
		for (int i = 0; i < n; i++) {
			if (w[i + n] > 0.0) {
				dx[i][k - 1] = 1;
			}
			else if (w[i + n] < 0.0) {
				dx[i][k - 1] = -1;
			}
			else {
				dx[i][k - 1] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		//% minimum
		MonotReduce(-1, A, dx[i], Ai);
		MonotReduce(-1, b, dx[i], bi);
		f(Ai, bi, xi);
		//% maximum
		MonotReduce(1, A, dx[i], Ai);
		MonotReduce(1, b, dx[i], bi);
		f(Ai, bi, xii);
		ww[i] = interval(xi[i].inf(), xii[i].sup());
	}
	return Success;
} //% >>> Monotonicity approach <<<

int
//%----------------------------------------------------------------------------
//% >>> Monotonicity approach <<<
//% based on augmented system matrix.
//% Gives the best results among all MAs, but is expensive.
//%----------------------------------------------------------------------------
MonotAugmentedI(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), //% function used to solve PILS
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //%revised affine right-hand vector
	ivector& ww) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), precond = 0, err;
	int K = AAFR::getlast(); //% number of initial parameters e_1, ..., e_K
	int niter = 0;
	aafrvector xa(n), x0(2 * n), xf(n), xff(n), bi(n), bg(2 * n);
	aafrmatrix Ai(n, n), Ag(2 * n, 2 * n);
	dvector pb(n);
	ivector w(2 * n);
	ovector dx(n);
	dmatrix pA(n, n);

	ww = ivector(n);
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	Ag.fill_in(0.0);
	bg.fill_in(0.0);
	for (int i = 0; i < n; i++) {
		bg[i] = b[i];
		for (int j = 0; j < n; j++) {
			Ag(i, j) = A(i, j);
			Ag(i + n, j + n) = A(i, j);
		}
	}
	//% computing derivatives over p_k
	for (int k = 1; k < K + 1; k++) {
		//% computing right-hand vector of system with derivatives			  
		PartAAFR(A, pA, k);
		PartAAFR(b, pb, k);
		for (int i = 0; i < n; i++) {
			bg[i + n] = pb[i];
			for (int j = 0; j < n; j++) {
				Ag(i + n, j) = pA(i, j);
			}
		}
		err = f(Ag, bg, x0, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		w = reduce(x0);
		for (int i = 0; i < n; i++) {
			if (w[i + n] > 0.0) {
				dx[i][k - 1] = 1;
			}
			else if (w[i + n] < 0.0) {
				dx[i][k - 1] = -1;
			}
			else {
				dx[i][k - 1] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		//% minimum
		MonotReduce(-1, A, dx[i], Ai);
		MonotReduce(-1, b, dx[i], bi);
		f(Ai, bi, xf, niter);
		//% maximum
		MonotReduce(1, A, dx[i], Ai);
		MonotReduce(1, b, dx[i], bi);
		err = f(Ai, bi, xff, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		ww[i] = interval(xf[i].reduce().inf(), xff[i].reduce().sup());
	}
	return Success;
} //% >>> Monotonicity approach <<<

int
//%----------------------------------------------------------------------------
//% >>> Monotonicity approach <<<
//% based on p-solution. Gives better results than standard MA.
//%----------------------------------------------------------------------------
MonotPsolution(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), //% method used to solve PILS
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& ww) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), precond = 0, err;
	int K = AAFR::getlast(); //% number of initial parameters e_1, ..., e_K
	int niter = 0;
	aafrvector xa(n), bb(n), x0(n), xf(n), xff(n), bi(n);
	dvector pb(n);
	ivector w(n), wb(n), wwin(n);
	ovector dx(n);
	aafrmatrix Ai(n, n);
	dmatrix pA(n, n);

	ww = ivector(n);
	err = f(A, b, xa, niter); //% xa = x(p) p-solution to the system
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	time_t start = clock();
	//% computing derivatives over p_k
	for (int k = 1; k < K + 1; k++) {
		//% computing right-hand vector of system with derivatives
		PartAAFR(A, pA, k);
		PartAAFR(b, pb, k);
		bb = pb - pA * xa;
		err = f(A, bb, x0, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		w = reduce(x0);
		for (int i = 0; i < n; i++) {
			if (w[i] > 0.0) {
				dx[i][k - 1] = 1;
			}
			else if (w[i] < 0.0) {
				dx[i][k - 1] = -1;
			}
			else {
				dx[i][k - 1] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		//% minimum
		MonotReduce(-1, A, dx[i], Ai);
		MonotReduce(-1, b, dx[i], bi);
		err = f(Ai, bi, xf, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		//% maximum
		MonotReduce(1, A, dx[i], Ai);
		MonotReduce(1, b, dx[i], bi);
		err = f(Ai, bi, xff, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		ww[i] = interval(xf[i].reduce().inf(), xff[i].reduce().sup());
		wwin[i] = interval(xf[i].reduce().sup(), xff[i].reduce().inf());
	}
	return Success;
} //% >>> Motonocity approach <<<

int
//%----------------------------------------------------------------------------
//% >>> Motonocity approach <<<
//% based on p-solution. In order to decrease the
//% computational time, all derivatives are computed using the method for
//% solving systems with multiple right-hand side. 
//%----------------------------------------------------------------------------
MonotPsolutionMRHS(int(*f)(const aafrmatrix&, const aafrvector&, aafrvector&, int&), //% function used to solve PILS
	int(*ff)(const aafrmatrix&, const aafrmatrix&, aafrmatrix&, int&),
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& ww) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = b.size(), precond = 0, err;
	int K = AAFR::getlast(); //% number of initial parameters e_1, ..., e_K
	int niter = 0;
	aafrvector xa(n), xf(n), xff(n), bi(n), bbi(n);
	aafrmatrix bb(n, K), x0(n, K);
	dvector pb(n);
	imatrix w(n, K);
	ovector dx(n);
	aafrmatrix Ai(n, n);
	dmatrix pA(n, n);

	ww = ivector(n);
	err = f(A, b, xa, niter); //% xa = x(p) p-solution to the system
	if (err != Success && err != NonVerifiedSolution) { return err; }
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	//% computing derivatives over p_k
	//% computing right-hand vector of system with derivatives
	for (int k = 1; k < K + 1; k++) {
		PartAAFR(A, pA, k);
		PartAAFR(b, pb, k);
		bbi = pb - pA * xa;
		for (int i = 0; i < n; i++) {
			bb(i, k - 1) = bbi[i];
		}
	}
	//% computing bounds for derivatives
	time_t start = clock();
	err = ff(A, bb, x0, niter);
	if (err != Success && err != NonVerifiedSolution) { return err; }
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < K; j++) {
			w(i, j) = x0(i, j).reduce();
		}
	}
	for (int i = 0; i < n; i++) {
		for (int k = 1; k < K + 1; k++) {
			if (w(i, k - 1) > 0.0) {
				dx[i][k - 1] = 1;
			}
			else if (w(i, k - 1) < 0.0) {
				dx[i][k - 1] = -1;
			}
			else {
				dx[i][k - 1] = 0;
			}
		}
	}
	//% time consuming loop
	//% it would be nice if we can decrease 
	//% the computational time of this loop
	for (int i = 0; i < n; i++) {
		//% minimum
		MonotReduce(-1, A, dx[i], Ai);
		MonotReduce(-1, b, dx[i], bi);
		err = f(Ai, bi, xf, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		//% maximum
		MonotReduce(1, A, dx[i], Ai);
		MonotReduce(1, b, dx[i], bi);
		err = f(Ai, bi, xff, niter);
		if (err != Success && err != NonVerifiedSolution) { return err; }
		ww[i] = interval(xf[i].reduce().inf(), xff[i].reduce().sup());
	}
	return Success;
} //% >>> Motonocity approach <<<

//%============================================================================
//% The methods below eventually must be replaced by guaranteed method
//% and then removed.
//%============================================================================

void
//%----------------------------------------------------------------------------
//% Motonocity approach based on p-solution
//% used in Global Optimization 
//% Remark: the method must be replaced in GO by guaranteed method
//%----------------------------------------------------------------------------
MonotPsolution(int ii,
	const ivector& p,  //% interval vector (parameters)
	const omvector& M, //% vector of real matrices
	const ovector& b,  //% vector of real vectors
	ivector& pmin,
	ivector& pmax)
{
	int n = b[0].size();
	int K = p.size();
	int niter = 0;
	omvector A(n + K + 1), B(K);
	ovector bb(n + K + 1);
	ovector dx(n);
	ivector w(n), wmin(n), wmax(n), px(K + n);
	dvector t1(n); //% center od [a]
	dvector t3(n); //% radius of [a]
	dmatrix t2(n, K); //% L matrix
	dvector x0(n);

	for (int k = 0; k < K; k++) {
		px[k] = p[k];
	}
	Rump(p, M, b, w, niter); //% enclosure for the solution of initial system
	for (int k = K; k < n + K; k++) {
		px[k] = w[k - K];
	}
	for (int k = 0; k < K + n + 1; k++) {
		A[k] = dmatrix(n, n);
		bb[k] = dvector(n);
	}
	//% creating matrix A
	for (int k = 0; k < K + 1; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[k](i, j) = M[k](i, j);
			}
		}
	}
	for (int k = K + 1; k < n + K + 1; k++) {
		A[k].fill_in(0.0); //% OK
	}
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	for (int k = 0; k < K; k++) {
		//%--------------------------------------------------------
		//% Main loop
		//% computing derivatives over p_k
		//%--------------------------------------------------------
		bb[0] = b[k + 1];
		for (int kk = 1; kk < K + 1; kk++) {
			bb[kk].fill_in(0.0);
		}
		for (int kk = K + 1; kk < n + K + 1; kk++) {
			for (int i = 0; i < n; i++) {
				bb[kk][i] = -M[k + 1](i, kk - K - 1);
			}
		}
		Rump(px, A, bb, w, niter);
		for (int i = 0; i < n; i++) {
			if (w[i] > 0.0) {
				dx[i][k] = 1;
			}
			else if (w[i] < 0.0) {
				dx[i][k] = -1;
			}
			else {
				dx[i][k] = 0;
			}
		}
	}
	for (int k = 0; k < K; k++) {
		if (dx[ii][k] == 1) {
			pmin[k] = p[k].inf();
			pmax[k] = p[k].sup();
			continue;
		}
		if (dx[ii][k] == -1) {
			pmin[k] = p[k].sup();
			pmax[k] = p[k].inf();
			continue;
		}
		if (dx[ii][k] == 0) {
			pmin[k] = p[k];
			pmax[k] = p[k];
			continue;
		}
	}
}

void
//%----------------------------------------------------------------------------
//% >>> Monotnicity approach <<<
//% The augmented system with x and dx/dp_k as unknowns is
//% solved using some method for solving PILS. The system has the form
//% | A(p)    0 | = | b(p) |
//% |d_A(p) A(p)| = |d_b(p)|
//%----------------------------------------------------------------------------
Monot(const ivector&  p, //% interval vector
	const omvector& M, //% vector of real matrices
	const ovector&  b, //% vector of interval vectors
	ivector&        ww)	//% interval solution vector
{
	int n = b[0].size(), m = 2 * n;
	int K = p.size();
	int niter = 0;
	dvector t1(m), t3(m), x0(m);
	omvector A(K + 1);
	ovector c(K + 1), dx(n);
	ivector w(m), wmin(n), wmax(n);
	dmatrix t2(m, K);

	ww = dvector(n);
	for (int k = 0; k < K + 1; k++) {
		A[k] = dmatrix(2 * n, 2 * n);
		c[k] = dvector(2 * n);
	}
	for (int k = 0; k < K + 1; k++) {
		//% computing derivatives over p_k
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[k](i, j) = M[k](i, j);
			}
			for (int j = n; j < m; j++) {
				A[k](i, j) = 0.0;
			}
			c[k][i] = b[k][i];
		}
		for (int i = n; i < m; i++) {
			for (int j = 0; j < m; j++) {
				A[k](i, j) = 0.0;
			}
			for (int j = n; j < m; j++) {
				A[k](i, j) = M[k](i - n, j - n);
			}
			c[k][i] = 0.0; //% this line fixed the error (the empty interval in the result)
		}
	}
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	for (int k = 0; k < K; k++) {
		//% computing derivatives over p_k
		for (int i = n; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[0](i, j) = M[k + 1](i - n, j);
			}
			c[0][i] = b[k + 1][i - n];
		}
		bool conv = kolevOuter(A, c, p, w, t1, t3, t2, x0, niter);
		for (int i = n; i < m; i++) {
			if (w[i] > 0.0) {
				dx[i - n][k] = 1;
			}
			else if (w[i] < 0.0) {
				dx[i - n][k] = -1;
			}
			else {
				dx[i - n][k] = 0;
			}
		}
	}
	ivector pmin(K), pmax(K);
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < K; k++) {
			if (dx[i][k] == 1) {
				pmin[k] = p[k].inf();
				pmax[k] = p[k].sup();
				continue;
			}
			if (dx[i][k] == -1) {
				pmin[k] = p[k].sup();
				pmax[k] = p[k].inf();
				continue;
			}
			if (dx[i][k] == 0) {
				pmin[k] = p[k];
				pmax[k] = p[k];
				continue;
			}
		}
		bool conv = kolevOuter(M, b, pmin, wmin, t1, t3, t2, x0, niter); //% MIN
		conv = kolevOuter(M, b, pmax, wmax, t1, t3, t2, x0, niter); //% MAX
		ww[i] = interval(wmin[i].inf(), wmax[i].sup());
	}
} //% >>> Monotonicity approach <<<

void
//%----------------------------------------------------------------------------
//% The motonocity approach for the PLP problem: 
//% the system is differentiated over p_k (k=1,...,K)
//% and the new system with derivatives as unknowns is obtained
//% Then, the augmented system with x and dx/dp_k as unknowns is
//% solved using some method for solving PILS
//% Finally, the derivative of l(x,p) over p_k is computed from
//% dx_i/dp_k and dc/dp_k (in this approach we assume that c=const
//% So, dc/dp_k=0
//%================================================================
//% Remark:
//% I think we cannot sum the derivatives if we don't have the hull form them
//% If we have only a crude estimation, then one derivative can be very large,
//% whereas the other not, and this one influences the overall result, which
//% can then be not correct
//%----------------------------------------------------------------------------
MonotPLP(int minmax,
	const dvector&	c,
	const ivector&	p,
	const omvector&	M,
	const ovector&	b,
	interval&		mw)
{
	int n = b[0].size(), m = 2 * n;
	int K = p.size();
	int niter = 0;
	dvector t1(m), t3(m), x0(m);
	ivector w(m), pp(K), ww(n);
	omvector A(K + 1);
	ovector bb(K + 1);
	ovector dx(n);
	dmatrix t2(m, K);

	for (int k = 0; k < K + 1; k++) {
		A[k] = dmatrix(2 * n, 2 * n);
		bb[k] = dvector(2 * n);
	}
	for (int k = 0; k < K + 1; k++) {
		//% compute derivatives over p_k
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[k](i, j) = M[k](i, j);
			}
			for (int j = n; j < m; j++) {
				A[k](i, j) = 0.0;
			}
			bb[k][i] = b[k][i];
		}
		for (int i = n; i < m; i++) {
			for (int j = 0; j < m; j++) {
				A[k](i, j) = 0.0;
			}
			for (int j = n; j < m; j++) {
				A[k](i, j) = M[k](i - n, j - n);
			}
			bb[k][i] = 0.0; //% this line fixed the error (the empty interval in the result)
		}
	}
	for (int i = 0; i < n; i++) {
		dx[i] = dvector(K);
	}
	for (int k = 0; k < K; k++) {
		//% compute derivatives over p_k
		for (int i = n; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[0](i, j) = M[k + 1](i - n, j);
			}
			bb[0][i] = b[k + 1][i - n];
		}
		bool conv = kolevOuter(A, bb, p, w, t1, t3, t2, x0, niter);
		//% here the derivatives dx/dp_k multiplied by c must be added
		interval dl;
		dl = 0.0;
		for (int i = n; i < m; i++) {
			dl = dl + c[i - n] * w[i];
		}
		if (minmax == -1) {
			if (dl > 0.0) {
				pp[k] = p[k].inf();
			}
			else if (dl < 0.0) {
				pp[k] = p[k].sup();
			}
			else {
				pp[k] = p[k];
			}
		}
		else {
			if (dl > 0.0) {
				pp[k] = p[k].sup();
			}
			else if (dl < 0.0) {
				pp[k] = p[k].inf();
			}
			else {
				pp[k] = p[k];
			}
		}
	}
	bool conv = kolevOuter(M, b, pp, ww, t1, t3, t2, x0, niter); //% minimum
	eIntPLP(minmax, c, x0, t1, t2, t3, mw);
}
