#include <conio.h>
#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../param_solver/mism.h"
#include "../param_solver/hladik.h"
#include "../param_solver/hansenbliekrohn.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/iterative.h"
#include "../param_solver/degrauwe.h"

#define EPSX		1.0e-6
#define MAXITER		100

bool
spectralRadius(
	int					op,	//% 0 - pre-cond, 1 - post-cond
	const aafrmatrix&	A,	//% revised affine matrix
	dmatrix&			B	//% radius of conditioned matrix A
)
//%---------------------------------------------------------------
//%
//% Computes and prints out matrix B=Ac^{-1}*A
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows();
	dmatrix		R(n, n);
	aafrmatrix	AA(n, n), Rf(n, n);
	aafrvector	b(n);

	mid(A, R); //% midpoint matrix
	std::cout << "A=[";
	if (inv(R)) {
		if (op == 0) {
			AA = R * A; //% pre-conditioning
		}
		else {
			AA = A * R; //% post-conditining
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = AA(i, j).reduce().rad();
				std::cout << B(i, j) << " ";
			}
			if (i < n - 1) {
				std::cout << ";" << std::endl;
			}
		}
		std::cout << "];" << std::endl << std::endl;
		std::cout << std::endl;
		return true;
	}
	return false;
}

void
iterInner(
	const aafrvector&	x,	  //% p-solution (revised affine vector)
	dvector&			innl, //% lower bound of IEH
	dvector&			innu  //% upper bound of IEH
)
//%--------------------------------------------------------------- 
//%
//% Computes inner solution 
//% from vector of reduced affine forms
//%
//%--------------------------------------------------------------- 
{
	int		n = x.size();
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
}

void
iterInner(
	const aafrmatrix&	x,	  //% p-solution (matrix of revised affine forms)
	dmatrix&			innl, //% lower bound if inner solution
	dmatrix&			innu  //% upper bound if inner solution
)
//%--------------------------------------------------------------- 
//%
//% Computes inner solution (multiple right-hand case)
//% from matrix of reduced affine forms 
//%
//%--------------------------------------------------------------- 
{
	int		n = x.num_rows();
	int		m = x.num_cols();
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
}

//%---------------------------------------------------------------
//%
//% Iterative methods for solving 
//% parametric interval linear systems
//%
//%---------------------------------------------------------------

bool
gsOuter(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised right-hand vector
	aafrvector&			w, //% p-solution (revised affine vector)
	int&				niter //% number of iterations
)
//%---------------------------------------------------------------
//%
//% Interval-affine Gauss-Seidel iteration
//% Iterative scheme: x=D^{-1}(b-(L+U)x)
//% Origin: Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% pre-conditioning with minpoint inverse
//%
//%--------------------------------------------------------------- 
//http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
{
	int			n = A.num_rows();
	aafrvector	w1(n), bb(n);
	cvector		xpidx(1);
	dvector		xinf(n), xsup(n), xpcfs(1);
	ivector		xp(n);
	aafrmatrix	AA(n, n);
	dmatrix		R(n, n);

	//% initial enclosure
	if (!ISM(A, b, xp, R)) return false; // O(n^3*K)
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { //% O(n)
		xpcfs[0] = xp[i].mid();
		w[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//% pre-conditioning
	AA = R * A; //% O(m*n^3*K) (m-cost of multiplication)
	bb = R * b; //% O(n^2*K)
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		w1 = w;
		//% entire loop O(K*log(K)*n^2), O(m*n^2) (m-cost of multiplication)
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) { // lower triangular
				if (AA(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
					s = s + AA(i, j) * w[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) { // upper triangular
				if (AA(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
					s = s + AA(i, j) * w1[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			w[i] = (bb[i] - s) / AA(i, i); // O(K*log(K))
		}
	} while (dist1(w, w1) > EPSX);
	return true;
}

bool
gsOuterII(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised right-hand vector
	aafrvector&			w, //% p-solution (revised affine vector)
	int&				niter //% number of iterations
)
//%---------------------------------------------------------------
//%
//% Interval-affine Gauss-Seidel iteration
//% Iterative scheme: x=D^{-1}(b-(L+U)x)
//% Origin: Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% post-conditioning with minpoint inverse
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows();
	aafrvector	w1(n), bb(n);
	cvector		xpidx(1);
	dvector		xinf(n), xsup(n), xpcfs(1);
	ivector		xp(n);
	aafrmatrix	AA(n, n);
	dmatrix		R(n, n), Am(n, n);

	//% initial enclosure of the post-conditioned system
	if (!ISMP0(A, b, xp, R)) return false; //% O(n^3*K)
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		w[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	AA = A * R; //% post-conditioning
	bb = b;
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		w1 = w;
		for (int i = 0; i < n; i++) { // all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) { // lower triangular
				if (AA(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
					s = s + AA(i, j) * w[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) { // upper triangular
				if (AA(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
					s = s + AA(i, j) * w1[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			w[i] = (bb[i] - s) / AA(i, i); // O(K*log(K))
		}
	} while (dist1(w, w1) > EPSX);
	w = R * w; //% final solution
	return true;
}

bool
gsOuter(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrmatrix&	B, //% revised affine right-hand matrix
	aafrmatrix&			x0, //% p-solution (revised affine matrix)
	int&				niter //% number of iterations
)
//%---------------------------------------------------------------
//%
//% Interval-affine Gauss-Seidel iteration with pre-conditioning
//% for multiple right-hand side
//%--------------------------------------------------------------- 
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% The inverse of the midpoint matrix is required here
//% to precondition the system
//% RAx=Rb, R=(A^c)^{-1}
//%
//%--------------------------------------------------------------- 
//%http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
{
	int			n = A.num_rows();
	int			m = B.num_cols();
	cvector		xpidx(1);
	dvector		xpcfs(1);
	aafrmatrix	x1(n, m), Bb(n, m), AA(n, n);
	dmatrix		xinf(n, m), xsup(n, m), R(n, n);
	imatrix		xp(n, m);

	//% initial enclosure obtained using ISM
	ISM(A, B, xp, R); // O(m*n^3*K) 
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		for (int j = 0; j < m; j++) {
			xpcfs[0] = xp(i, j).mid();
			x0(i, j) = AAFR(xpidx, xpcfs, xp(i, j).rad());
		}
	}
	//% preconditioning with the midpoint inverse
	AA = R * A; // O(n^3*K)
	Bb = R * B; // O(n^2*K)
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		imatrix xi(n, m);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				xi(i, j) = x1(i, j).reduce(); //% O(n)
			}
		}
		for (int i = 0; i < n; i++) { //% all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			for (int k = 0; k < m; k++) {
				AAFR s = 0.0;
				for (int j = 0; j < i; j++) {
					if (AA(i, j).reduce() != 0.0 && x0(j, k).reduce() != 0.0) {
						s = s + AA(i, j) * x0(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				for (int j = i + 1; j < n; j++) {
					if (AA(i, j).reduce() != 0.0 && x1(j, k).reduce() != 0.0) {
						s = s + AA(i, j) * x1(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
					}
				}
				x0(i, k) = (Bb(i, k) - s) / AA(i, i); //% O(K*log(K))
			}
		}
	} while (dist1(x0, x1) > EPSX);
	return true;
}

bool
Krawczyk(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised affine vector
	aafrvector&			w, //% p-solution (revised affine vector)
	int&				niter //% number of iterations
)
//%---------------------------------------------------------------
//%
//% Krawczyk iteration wihtout residual correction
//% the starting vector must be computed by the
//% method without residual correction 
//%
//%---------------------------------------------------------------
{
	int	n = A.num_rows();
	cvector idx(1);
	dvector	cfs(1);
	dmatrix	R(n, n), I(n, n);
	ivector	w0(n);
	aafrmatrix V(n, n), B(n, n);
	aafrvector v(n), y0(n), y1(n);

	idx[0] = 0;
	if (HansenBliekRohn(A, b, w0, R)) { //% initial solution
		V = R * A; //% preconditioning of A
		v = R * b; //% preconditioning of b
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - V(i, j); //% (I-V(e))
			}
		}
		for (int i = 0; i < n; i++) {
			cfs[0] = w0[i].mid();
			y0[i] = AAFR(idx, cfs, w0[i].rad()); //% starting point
		}
		//%-------------------------------------------------------------
		//% ITERATION
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
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
		} while (distH(y0, y1) > EPSX);
		w = y0;
		return true;
	}
	return false;
}

bool
Krawczyk(
	const aafrmatrix&	A,
	const aafrmatrix&	bb,
	aafrmatrix&			w,
	int&				niter
)
//%---------------------------------------------------------------
//%
//% Krawczyk iteration
//% from paper with Milan
//% Multiple right-hand side
//%
//%---------------------------------------------------------------
{
	int			n = A.num_rows();
	int			m = bb.num_cols();
	real		eps = 0.01, mi = 1.0e-15;
	cvector		idx(1);
	dvector		cfs(1);
	dmatrix		bm(n, m), x0(n, m), R(n, n), I(n, n);
	imatrix		w0(n, m);
	aafrmatrix	V(n, n), B(n, n), vv(n, m), y0(n, m), y1(n, m);

	idx[0] = 0;
	if (BauerSkeel2_MRHS(A, B, w0, R)) {
		//if (HansenBliekRohn(A, b, w0, R)) { // OK
		mid(bb, bm); // midpoint of b
		x0 = R * bm; // midpoint solution
		V = R * A; // preconditioning of A
		vv = R * bb; // preconditioning of b
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - V(i, j); // matrix (I-V(e))
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				cfs[0] = w0(i, j).mid();
				y0(i, j) = AAFR(idx, cfs, w0(i, j).rad()); // starting point
			}
		}
		//%-------------------------------------------------------------
		//% ITERATION
		//%-------------------------------------------------------------
		niter = 0;
		do {
			niter++;
			if (niter > MAXITER) break;
			y1 = y0;
			for (int i = 0; i < n; i++) { //% all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
				for (int k = 0; k < m; k++) {
					AAFR s = 0.0;
					for (int j = 0; j < n; j++) {
						if (B(i, j).reduce() != 0.0 && y0(j, k).reduce() != 0.0) {
							s = s + B(i, j) * y0(j, k); //% O(K*log(K)), O(m), m-cost of multiplication
						}
					}
					y0(i, k) = vv(i, k) + s; //% O(K*log(K))
				}
			}

		} while (dist1(y0, y1) > EPSX);
		w = y0;
		return true;
	}
	return false;
}

bool
mkolevRC(
	const aafrmatrix&	A,	//% in: system matrix
	const aafrvector&	b,	//% in: right-hand vector
	aafrvector&			w,	//% out: solution
	int&				niter //% out: number of iterations
)
//%--------------------------------------------------------------- 
//%
//% Modified Kolev's Iteration method 
//% with residual correction
//% based on revised affine forms
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows(); //% problem size
	aafrvector	d(n), v(n), v1(n);
	cvector		xpidx(1);
	dvector		b0(n), x0(n), xpcfs(1);
	ivector		v2(n), xp(n);
	aafrmatrix	C(n, n), C0(n, n);
	dmatrix		A0(n, n), R(n, n);

	mid(b, b0); //% midpoint vector mid(b)
	ISM(A, b, xp, R); //% xp - initial enclosure
	x0 = R * b0; //% midpoint solution
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C0(i, j) = A(i, j) - A(i, j).cfs()[0];
		}
	}
	d = R * (b - A * x0); //% residual correction
	C = R * C0;
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { //% O(n)
		xpcfs[0] = xp[i].mid();
		v[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//%-------------------------------------------------------------
	//% ITERATION
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
	for (int i = 0; i < n; i++) {
		w[i] = (x0[i] + v[i]);
	}
	return true;
}

bool
mkolevRCII(
	const aafrmatrix&	A, //% in: system matrix
	const aafrvector&	b,	//% in: right-hand vector
	aafrvector&			w,	//% out: solution
	int&				niter //% out: number of iterations
)
//%--------------------------------------------------------------- 
//%
//% Modified Kolev's Iteration method 
//% with residual correction and post-conditioning
//% based on revised affine forms
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows(); //% problem size
	aafrvector	d(n), v(n), v1(n);
	cvector		xpidx(1);
	dvector		b0(n), x0(n), xpcfs(1);
	ivector		xp(n);
	aafrmatrix	C(n, n), C0(n, n);
	dmatrix		A0(n, n), R(n, n);

	mid(b, b0); //% midpoint vector mid(b)
	ISMP0(A, b, xp, R); //% xp - initial enclosure
	x0 = R * b0; //% midpoint solution
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C0(i, j) = A(i, j) - A(i, j).cfs()[0];
		}
	}
	C = C0 * R; //% post-conditioning
	d = (b - A * x0); //% residual correction
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { //% O(n)
		xpcfs[0] = xp[i].mid();
		v[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//%-------------------------------------------------------------
	//% ITERATION
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
	return true;
}

bool
mkolev(
	const aafrmatrix&	A,	//% in: system matrix
	const aafrvector&	b,	//% in: right-hand vector
	aafrvector&			w,	//% out: solution
	int& niter				//% out: number of iterations
)
//%--------------------------------------------------------------- 
//%
//% Modified Kolev's Iteration method 
//% without residual correction
//% based on revised affine forms
//%
//%--------------------------------------------------------------- 
{
	int				n = A.num_rows();
	aafrmatrix		C(n, n), C0(n, n);
	aafrvector		d(n), bb(n), x0f(n), v(n), v1(n);
	aafrvector		dd(n);
	dmatrix			A0(n, n), R(n, n);
	dvector			b0(n), x0(n);
	ivector			v2(n), xp(n);
	cvector			xpidx(1);
	dvector			xpcfs(1);

	ISM(A, b, xp, R); //% xp - initial enclosure
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C0(i, j) = A(i, j) - A(i, j).cfs()[0];
		}
	}
	d = R * b; //% no residual correction
	C = R * C0; //% pre-conditioning
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		v[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		v1 = v;
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				if (C(i, j).reduce() != 0.0 && v[j].reduce() != 0.0) {
					s = s + C(i, j) * v[j];
				}
			}
			v[i] = d[i] - s; //% update entry
		}
	} while (dist1(v1, v) > EPSX);
	w = v;
	return true;
}

bool
mkolevII(
	const aafrmatrix&	A,
	const aafrvector&	b,
	aafrvector&			w,
	int& niter
)
//%--------------------------------------------------------------- 
//% Modified version of Kolev's method
//% computing with reduced affine forms
//% with post-conditioning
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows(), niter2;
	aafrvector	d(n), x0f(n), v(n), v1(n);
	cvector		xpidx(1);
	dvector		x0(n), xpcfs(1);
	ivector		xp(n);
	aafrmatrix	C(n, n), AA(n, n), Af(n, n);
	dmatrix		R(n, n);

	ISMP0(A, b, xp, R); //% initial enclosure
	AA = A * R;
	d = b;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C(i, j) = AA(i, j) - AA(i, j).cfs()[0];
		}
	}
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		v[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
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
	} while (dist1(v, v1) > EPSX);
	for (int i = 0; i < n; i++) {
		w[i] = v[i];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Af(i, j) = A(i, j).cfs()[0];
		}
	}
	mkolev(Af, v, w, niter2);
	return true;
}

//%===============================================================

bool
jacobi(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised affine right-hand vector
	aafrvector&			x0, //% p-solution (revised affine vector)
	int&				niter //% number of iterations
)
//%---------------------------------------------------------------
//%
//% Jacobi iteration
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% RAx=Rb, R=(A^c)^{-1}
//%
//%--------------------------------------------------------------- 
{

	int			n = A.num_rows();
	aafrvector	x1(n), D(n), bb(n);
	cvector		xpidx(1);
	dvector		xinf(n), xsup(n), xpcfs(1);
	ivector		xp(n);
	aafrmatrix	AA(n, n), Rf(n, n);
	dmatrix		R(n, n);

	//% Initial enclosure computed using ISM
	ISM(A, b, xp, R);
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) {
		xpcfs[0] = xp[i].mid();
		x0[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Rf(i, j) = AAFR(R(i, j));
		}
	}
	AA = Rf * A;
	bb = Rf * b;
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				if (i != j) {
					if (AA(i, j).reduce() != 0.0 && x1[j].reduce() != 0.0) {
						s = s + AA(i, j) * x1[j];
					}
				}
			}
			x0[i] = (bb[i] - s) / AA(i, i);
		}
	} while (dist1(x0, x1) > EPSX);
	return true;
}

bool
gsOuter0(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised right-hand vector
	aafrvector&			x0, //% interval solution
	int&				niter //% number of iterations
)
//%--------------------------------------------------------------- 
//%
//% Parametric Gauss-Seidel iteration
//% Iterative scheme: x=D^{-1}(b-(L+U)x)
//% Origin: Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% (preconditioning with minpoint inverse)
//% the difference is that 
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows();
	aafrvector	bb(n), x1(n);
	cvector		xpidx(1);
	dvector		xinf(n), xsup(n), xpcfs(1);
	ivector		xp(n), xi(n), xii(n);
	aafrmatrix	AA(n, n);
	dmatrix		R(n, n);

	//% initial enclosure obtained using ISM
	ISM(A, b, xp, R); // O(n^3*K)
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		x0[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	//% preconditioning with the midpoint inverse
	AA = R * A; // O(n^3*K)
	bb = R * b; // O(n^2*K)
	//%-------------------------------------------------------------
	//% ITERATION
	//%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
		x1 = x0;
		xi = reduce(x1); // O(n)
		for (int i = 0; i < n; i++) { // all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) {
				if (AA(i, j).reduce() != 0.0 && x0[j].reduce() != 0.0) {
					s = s + AA(i, j) * x0[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) {
				if (AA(i, j).reduce() != 0.0 && x1[j].reduce() != 0.0) {
					s = s + AA(i, j) * x1[j]; // O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			xii[i] = ((bb[i] - s) / AA(i, i)).reduce();
			xii[i] = interval(ISLIB_MAX(xi[i].inf(), xii[i].inf()), ISLIB_MIN(xi[i].sup(), xii[i].sup()));
			xpcfs[0] = xii[i].mid();
			x0[i] = AAFR(xpidx, xpcfs, xii[i].rad());
		}
	} while (distH(xi, xii) > EPSX);
	return true;
}

bool
sor(
	const aafrmatrix&	A, //% revised affine matrix
	const aafrvector&	b, //% revised right-hand vector
	aafrvector&			x0, //% interval solution
	real				wgt, //% relaxation parameter
	int&				niter //% number of iterations
)
//%--------------------------------------------------------------- 
//%
//% Succesive overrelaxation method
//% Iterative method (preconditioning with minpoint inverse)
//%
//%--------------------------------------------------------------- 
{
	int			n = A.num_rows();
	cvector		xpidx(1);
	dvector		xpcfs(1);
	ivector		xp(n);
	dmatrix		R(n, n);
	aafrvector	x1(n), bb(n);
	aafrmatrix	AA(n, n), Rf(n, n);

	ISM(A, b, xp, R); //% initial enclosure
	//% transform initial interval enclosure into revised affine form
	xpidx[0] = 0;
	for (int i = 0; i < n; i++) { // O(n)
		xpcfs[0] = xp[i].mid();
		x0[i] = AAFR(xpidx, xpcfs, xp[i].rad());
	}
	AA = R * A;
	bb = R * b;
	//%-------------------------------------------------------------
	//% ITERATION
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
	return true;
}

bool
jacobi(
	const imatrix&	A,
	const ivector&	b,
	ivector&		w
)
//%--------------------------------------------------------------- 
//%
//% Jacobi iteration for standard interval linear systems
//% without taking into account dependencies 
//% TRZEBA PRZEROBIC NA WERSJE "ELEMENTOWA"
//%
//%---------------------------------------------------------------
{
	real EPS = 1.0e-9, error;
	int n = A.num_rows();
	dvector xinf(n), xsup(n);
	dmatrix R(n, n);
	ivector x0(n), x1(n);
	ivector D(n), bb(n);
	imatrix L(n, n), U(n, n), Dinv(n, n), AA(n, n);

	x0.fill_in(0.0); // initial point
	L.fill_in(0.0);
	U.fill_in(0.0);
	Dinv.fill_in(0.0);
	mid(A, R);
	inv(R);
	AA = R * A;
	bb = R * b;
	for (int i = 0; i < n; i++) {
		D[i] = AA(i, i);
		Dinv(i, i) = 1.0 / D[i];
		for (int j = 0; j < i; j++) {
			L(i, j) = AA(i, j);
			U(j, i) = AA(j, i);
		}
	}
	do {
		error = 0.0;
		x1 = x0;
		x0 = Dinv * (bb - (L + U)*x0);
		for (int i = 0; i < n; i++) {
			error += ISLIB_ABS(x0[i].inf() - x1[i].inf()) + ISLIB_ABS(x0[i].sup() - x1[i].sup());
		}
	} while (error > EPS);
	w = x0;
	return true;
}

bool
RohnOverdeterm(
	const aafrmatrix& A,
	const aafrvector& b,
	ivector& w
)
{
	int n = A.num_rows(), m = A.num_cols();
	dmatrix Ac(n, m);
	mid(A, Ac);
	aafrmatrix Ap = Ac * A;
	aafrvector bp = Ac * b;

	ISM(Ap, bp, w);
}