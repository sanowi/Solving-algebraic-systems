#include <iostream>
#include <conio.h>
#include <omp.h>
#include <ctime>

#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/gausselim.h"
#include "../utils/inverse.h"
#include "../newton/Newton.h"
#include "../affine/aaf.h"
#include "../affine/aafr.h"
#include "../utils/randnum.h"
#include "kolev.h"

#define EPSP 1.0e-16
#define MAX_ITER 500

void
kolevReducei(int i,
	const dvector& t1,
	const dmatrix& t2,
	const dvector& t3,
	interval& v)
	//%----------------------------------------------------------------------------
	//% reduces i-th element of the parametric solution to an interval
	//%----------------------------------------------------------------------------
{
	interval p(-1.0, 1.0);
	int n = t2.num_rows();
	int K = t2.num_cols();
	interval s(0.0);
	for (int k = 0; k < K; k++) {
		s = s + interval(t2(i, k)) * p; // t2 <-> L
	}
	v = interval(t1[i]) + s + interval(t3[i]) * interval(-1.0, 1.0); //% t1 + (s + t3)[-1,1]
}

void
//%----------------------------------------------------------------------------
//% reduce the parametric solution to an interval (quadratic)
//%----------------------------------------------------------------------------
kolevReduceQ(const dvector& t1,
	const dmatrix& t2,
	const dvector& t3,
	const dmatrix& t4,
	ivector& v)
{
	interval p(-1.0, 1.0);
	interval p2(0.0, 1.0);
	int n = t2.num_rows();
	int K = t2.num_cols();
	for (int i = 0; i < n; i++) {
		interval s(0.0);
		interval s2(0.0);
		for (int k = 0; k < K; k++) {
			s = s + t2(i, k) * p; // t2 = L in Kolev's paper: L*p
			s2 = s2 + t4(i, k) * p2; // t4 = Q in Kolev's paper: Q*p^2
		}
		v[i] = t1[i] + s + s2 + t3[i] * p; // t1 + (s + t3)[-1,1]
	}
}

void
//%----------------------------------------------------------------------------
//% Transforms the form M1=A0+sum(Ai*pi), b1=b0+sum(bi*pi)
//% to the form with omatrix oA and ovector ob used for trusses
//%----------------------------------------------------------------------------
to_oA_ob(const omvector& M1,
	const ovector& b1,
	int K,
	omatrix& oA,
	ovector& ob)
{
	int i, j, k, n = b1.size();
	dvector el(K);

	oA = omatrix(n, n);
	ob = ovector(n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				el[k] = M1[k](i, j);
			}
			oA(i, j) = el;
		}
		for (k = 0; k < K; k++) {
			el[k] = b1[k][i];
		}
		ob[i] = el;
	}
}

void
//%----------------------------------------------------------------------------
//% Computes A(p) from A^{(k)} (stored in M) and p
//%----------------------------------------------------------------------------
omvec2real(const omvector& M,
	const dvector& p,
	dmatrix& A)
{
	int n = M[0].num_rows();
	int K = p.size();
	A = dmatrix(n, n); A = M[0];
	for (int k = 1; k < K + 1; k++) {
		A = A + M[k] * p[k - 1];
	}
}

void
//%----------------------------------------------------------------------------
//% Computes b(p) from b^{(k)} (stored in b) and p
//%----------------------------------------------------------------------------
ovec2real(const ovector& b,
	const dvector& p,
	dvector& bb)
{
	int n = b[0].size();
	int K = p.size();
	bb = dvector(n); bb = b[0];
	for (int k = 1; k < K + 1; k++) {
		bb = bb + b[k] * p[k - 1];
	}
}

real
//%----------------------------------------------------------------------------
//% stopping criterion for Kolev's method:
//% q([x],[y])=max|xi-yi|
//%----------------------------------------------------------------------------
stopCrit(const ivector& x,
	const ivector& y)
{
	int n = x.size();
	real s = 0.0;
	for (int i = 0; i < n; i++) {
		s = fmax(s, fmax(fabs(x[i].inf() - y[i].inf()), fabs(x[i].sup() - y[i].sup())));
	}
	return s;
}

void
//%----------------------------------------------------------------------------
//% reduce the parametric solution to an interval solution (interval vector)
//%----------------------------------------------------------------------------
reduce(const dvector& t1,
	const dmatrix& t2,
	const dvector& t3,
	ivector& v)
{
	interval p(-1.0, 1.0);
	int n = t2.num_rows();
	int K = t2.num_cols();
	for (int i = 0; i < n; i++) {
		interval s(0.0);
		for (int k = 0; k < K; k++) {
			s = s + t2(i, k) * p; // t2 is the same as L in Kolev's paper (L*p)
		}
		v[i] = t1[i] + s + t3[i] * interval(-1.0, 1.0); // t1 + (s + t3)[-1,1]
	}
}

void
//%----------------------------------------------------------------------------
//% transforms the original system to the system with [p_i]=[-1,1]
//%----------------------------------------------------------------------------
ktransform(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	omvector& M,
	ovector& b,
	ivector& p)
{
	int n = b1[0].size();
	int K = p1.size();
	for (int i = 0; i < K + 1; i++) {
		M[i] = dmatrix(n, n);
		b[i] = dvector(n);
	}
	M[0] = M1[0];
	b[0] = b1[0];

	for (int k = 1; k < K + 1; k++) {
		M[0] = M[0] + M1[k] * p1[k - 1].mid();
		b[0] = b[0] + b1[k] * p1[k - 1].mid();
	}
	for (int k = 1; k < K + 1; k++) {
		M[k] = M1[k] * p1[k - 1].rad();
		b[k] = b1[k] * p1[k - 1].rad();
		p[k - 1] = interval(-1.0, 1.0);
	}
}

void
//%----------------------------------------------------------------------------
//% transforms the original system to the system with [p_i]=[-1,1]
//%----------------------------------------------------------------------------
ktransform(const omvector& M1,
	const ivector& p1,
	omvector& M,
	ivector& p)
{
	int n = M1[0].num_rows();
	int K = p1.size();
	for (int i = 0; i < K + 1; i++) {
		M[i] = dmatrix(n, n);
	}
	M[0] = M1[0];

	for (int k = 1; k < K + 1; k++) {
		M[0] = M[0] + M1[k] * p1[k - 1].mid();
	}
	for (int k = 1; k < K + 1; k++) {
		M[k] = M1[k] * p1[k - 1].rad();
		p[k - 1] = interval(-1.0, 1.0);
	}
}

void
//%----------------------------------------------------------------------------
//% transforms the original system to the system with [p_i]=[-1,1]
//%----------------------------------------------------------------------------
ktransformEig(const omvector& M1,
	const ivector& p1,
	omvector& M,
	ivector& p)
{
	int n = M1[0].num_rows();
	int K = p1.size();
	for (int k = 0; k < K + 1; k++) {
		M[k] = M1[k];
	}

	for (int k = 1; k < K + 1; k++) {
		M[0] = M[0] + M1[k] * p1[k - 1].mid();
	}
	for (int k = 0; k < K; k++) {
		real r = p1[k].rad();
		p[k] = interval(-r, r);
	}
}

void
//%---------------------------------------------------------------------------- 
//% calculates only inner solution
//% which is returned in variable innsol
//%----------------------------------------------------------------------------
kolevInner(const dvector& x0,
	const dvector& t1,
	const dvector& t3,
	const dmatrix& t2,
	ivector& innsol)
{
	int n = t1.size();
	int K = t2.num_cols();
	dmatrix mt2(n, K);
	ivector a(n);
	dvector au(n), al(n), pv(K), innu(n), innl(n);
	mag(t2, mt2);
	pv.fill_in(1.0);
	pv = mt2 * pv;
	a = x0 + t1 + t3 * interval(-1.0, 1.0);
	sup(a, au);
	inf(a, al);
	innl = au - pv;
	innu = al + pv;
	for (int i = 0; i < n; i++) {
		innsol[i] = interval(innl[i], innu[i]);
	}
}

void
//%----------------------------------------------------------------------------
//% Version of KolevInner
//% calculates the endpoints of the inner solution: innl (left) and innu (right)
//% there may be the situation that innl >= innu
//%----------------------------------------------------------------------------
kolevInner(const dvector& x0,
	const dvector& t1,
	const dvector& t3,
	const dmatrix& t2,
	dvector& innl,
	dvector& innu)
{
	int n = t1.size();
	int K = t2.num_cols();
	dmatrix mt2(n, K);
	ivector a(n);
	dvector au(n), al(n), pv(K);

	innl = dvector(n);
	innu = dvector(n);
	mag(t2, mt2);
	pv.fill_in(1.0);
	pv = mt2 * pv;
	a = x0 + t1 + t3 * interval(-1.0, 1.0);
	sup(a, au);
	inf(a, al);
	innl = au - pv;
	innu = al + pv;
}

bool
//%----------------------------------------------------------------------------
//% Original Kolev's method; calculates only outer solution
//% it is assumed here that the parameters vary within the interval [-1,1]
//%----------------------------------------------------------------------------
kolevOuter0(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w,
	dvector& t1,
	dvector& t3,
	dmatrix& t2,
	dvector& x0,
	int& niter)
{
	int n = b1[0].size();
	int K = p1.size();
	real eps = 1.0e-6;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n);

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	ktransform(M1, b1, p1, M, b, p);
	A0 = M[0]; // midpoint matrix of the new system
	b0 = b[0]; // midpoint vector of the new system
	R = A0;
	x0 = dvector(n);
	if (inv(R)) { // if the mipoint matrix is nonsingular, the computation continues
		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), D0(K);
		ovector d(K), G(K);

		//gausse(A0, b0, x0); // computing x0 - the midpoint solution
		x0 = R * b0;
		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); // d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; // B^k=R*A^k, tez chyba OK
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		// v^(0)=0
		v = C0 * p; // v^(1)=C*p 

		// v^(2)=t1+t2*p+t3	
		t1 = dvector(n);
		t3 = dvector(n);
		t2 = dmatrix(n, K);

		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);

		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0;
		}

		for (int i = 0; i < n; i++) {
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					s += fabs(D[k](i, j) + D[j](i, k));
				}
				t1[i] = t1[i] - D[k](i, k) / 2.0;
				t2(i, k) = C0(i, k); // t2 <-> L
				t3[i] = t3[i] + fabs(D[k](i, k) / 2.0);
			}
			t3[i] = t3[i] + s;
		}
		niter = 0;
		////////////////////////////////////////////////////////////////////////
		//% ITERATION: computing v^(k), k>=3
		//%
		reduce(t1, t2, t3, v);
		do {
			if (niter > MAX_ITER) return false;
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				G[k] = B[k] * t1;
				for (int j = 0; j < n; j++) {
					C(j, k) = C0(j, k) - G[k][j];
				}
				D0[k] = B[k] * t2;
			}
			t1.fill_in(0.0);
			t3i = t3;
			t3.fill_in(0.0);
			for (int i = 0; i < n; i++) { // n is the size of v
				real s = 0.0;
				real s1 = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < n; j++) {
						s1 += fabs(B[k](i, j)) * t3i[j];
					}
					for (int j = k + 1; j < K; j++) {
						s += fabs(D0[k](i, j) + D0[j](i, k));
					}
					t1[i] = t1[i] - D0[k](i, k) / 2.0;
					t2(i, k) = C(i, k); // this is the L operator, which is multiplied by p
					t3[i] = t3[i] + fabs(D0[k](i, k) / 2.0);
				}
				t3[i] = t3[i] + s + s1; // the p-solution is in the form: t1 + t2*p + [t3]
			}
			reduce(t1, t2, t3, v);
		} while (stopCrit(v, v1) > eps);
		reduce(t1, t2, t3, v);
		for (int i = 0; i < n; i++) { // final result for outer solution
			w[i] = v[i] + interval(x0[i]);
		}
		//std::cout << "Kolev's method: Iterations: " << niter << std::endl;
		return true;
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Original Kolev's method; calculates only outer solution
//% it is assumed here that the parameters vary within the interval [-1,1]
//%----------------------------------------------------------------------------
kolevOuter(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w,
	dvector& t1,
	dvector& t3,
	dmatrix& t2,
	dvector& x0,
	int& niter)
{
	int n = b1[0].size();
	int K = p1.size();
	real eps = 1.0e-16;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n);
	imatrix RI(n, n);

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	x0 = dvector(n);
	t1 = dvector(n);
	t3 = dvector(n);
	t2 = dmatrix(n, K);

	t1.fill_in(0.0);
	t2.fill_in(0.0);
	t3.fill_in(0.0);

	ktransform(M1, b1, p1, M, b, p);
	A0 = M[0]; //% midpoint matrix of the new system
	b0 = b[0]; //% midpoint vector of the new system
	R = A0;
	to_int(R, RI);
	x0 = dvector(n);
	if (inv(R)) { //% if the mipoint matrix is nonsingular, the computation continues
		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), D0(K);
		ovector d(K), G(K);
		gausse(A0, b0, x0); //% computing x0 - the midpoint solution
		//% x0 = R * b0;
		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); //% d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; //% B^k=R*A^k, tez chyba OK
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		//% v^(0)=0
		v = C0 * p; //% v^(1)
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0;
		}
		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);

		for (int i = 0; i < n; i++) {
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					s += fabs(D[k](i, j) + D[j](i, k));
				}
				t1[i] = t1[i] - D[k](i, k) / 2.0;
				t2(i, k) = C0(i, k); // t2 <-> L
				t3[i] = t3[i] + fabs(D[k](i, k) / 2.0);
			}
			t3[i] = t3[i] + s;
		}
		int niter = 0;
		////////////////////////////////////////////////////////////////////////
		//% ITERATION: computing v^(k), k>=3
		//%
		reduce(t1, t2, t3, v);
		do {
			if (niter > MAX_ITER) {
				//std::cout << "Exited with false" << std::endl;
				return false;
			}
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				G[k] = B[k] * t1;
				for (int j = 0; j < n; j++) {
					C(j, k) = C0(j, k) - G[k][j];
				}
				D0[k] = B[k] * t2;
			}
			t1.fill_in(0.0);
			t3i = t3;
			t3.fill_in(0.0);
			for (int i = 0; i < n; i++) { //% n is the size of v
				real s = 0.0;
				real s1 = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < n; j++) {
						s1 += fabs(B[k](i, j)) * t3i[j];
					}
					for (int j = k + 1; j < K; j++) {
						s += fabs(D0[k](i, j) + D0[j](i, k));
					}
					t1[i] = t1[i] - D0[k](i, k) / 2.0;
					t2(i, k) = C(i, k); //% this is the L operator, which is multiplied by p
					t3[i] = t3[i] + fabs(D0[k](i, k) / 2.0);
				}
				t3[i] = t3[i] + s + s1; //% the p-solution is in the form: t1 + t2*p + [t3]
			}
			reduce(t1, t2, t3, v);
		} while (stopCrit(v, v1) > EPSP);
		reduce(t1, t2, t3, v);
		for (int i = 0; i < n; i++) { //% final result for outer solution
			w[i] = v[i] + interval(nextafter(x0[i], -ISLIB_INFINITY), nextafter(x0[i], ISLIB_INFINITY));
		}
		return true;
	}
	std::cout << "Exited with false" << std::endl;
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Original Kolev's method; calculates only outer solution
//% it is assumed here that the parameters vary within the interval [-1,1]
//%----------------------------------------------------------------------------
kolevOuter(int mode,
	const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w)
{
	int n = b1[0].size();
	int K = p1.size();
	real eps = 1.0e-10; //% 1.0e-16;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n);
	imatrix RI(n, n);

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	dvector t1(n), t3(n), x0(n);
	dmatrix t2(n, K);

	if (mode == 0) {
		RoundDown();
	}
	else {
		RoundUp();
	}
	ktransform(M1, b1, p1, M, b, p);
	A0 = M[0]; //% midpoint matrix of the new system
	b0 = b[0]; //% midpoint vector of the new system
	R = A0;
	to_int(R, RI);
	x0 = dvector(n);
	if (inv(R)) { //% if the mipoint matrix is nonsingular, the computation continues
		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), D0(K);
		ovector d(K), G(K);
		gausse(A0, b0, x0); //% computing x0 - the midpoint solution
		//% x0 = R * b0;
		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); //% d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; //% B^k=R*A^k, tez chyba OK
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		//% v^(0)=0
		v = C0 * p; //% v^(1)
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0;
		}
		//% v^(2)=t1+t2*p+t3	
		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);

		for (int i = 0; i < n; i++) {
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					s += fabs(D[k](i, j) + D[j](i, k));
				}
				t1[i] = t1[i] - D[k](i, k) / 2.0;
				t2(i, k) = C0(i, k); // t2 <-> L
				t3[i] = t3[i] + fabs(D[k](i, k) / 2.0);
			}
			t3[i] = t3[i] + s;
		}
		int niter = 0;
		////////////////////////////////////////////////////////////////////////
		//% ITERATION: computing v^(k), k>=3
		//%
		reduce(t1, t2, t3, v);
		do {
			if (niter > MAX_ITER) {
				std::cout << "Exited with false" << std::endl;
				return false;
			}
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				G[k] = B[k] * t1;
				for (int j = 0; j < n; j++) {
					C(j, k) = C0(j, k) - G[k][j];
				}
				D0[k] = B[k] * t2;
			}
			t1.fill_in(0.0);
			t3i = t3;
			t3.fill_in(0.0);
			for (int i = 0; i < n; i++) { //% n is the size of v
				real s = 0.0;
				real s1 = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < n; j++) {
						s1 += fabs(B[k](i, j)) * t3i[j];
					}
					for (int j = k + 1; j < K; j++) {
						s += fabs(D0[k](i, j) + D0[j](i, k));
					}
					t1[i] = t1[i] - D0[k](i, k) / 2.0;
					t2(i, k) = C(i, k); //% this is the L operator, which is multiplied by p
					t3[i] = t3[i] + fabs(D0[k](i, k) / 2.0);
				}
				t3[i] = t3[i] + s + s1; //% the p-solution is in the form: t1 + t2*p + [t3]
			}
			reduce(t1, t2, t3, v);
		} while (stopCrit(v, v1) > EPSP);
		reduce(t1, t2, t3, v);
		for (int i = 0; i < n; i++) { //% final result for outer solution
			w[i] = v[i] + interval(nextafter(x0[i], -ISLIB_INFINITY), nextafter(x0[i], ISLIB_INFINITY));
		}
		std::cout << "Exited with true" << std::endl;
		return true;
	}
	std::cout << "Exited with false" << std::endl;
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Original Kolev's method; calculates only outer solution
//% it is assumed here that the parameters vary within the interval [-1,1]
//%----------------------------------------------------------------------------
kolevOuter(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w)
{
	int n = M1[0].num_rows();
	ivector w1(n), w2(n);
	bool conv = kolevOuter(0, M1, b1, p1, w1);
	if (!conv) {
		RoundNear();
		return false;
	}
	//std::cout << w1 << std::endl;
	conv = kolevOuter(1, M1, b1, p1, w2);
	if (!conv) {
		RoundNear();
		return false;
	}
	//std::cout << w2 << std::endl;
	RoundNear();
	w = ivector(n);
	for (int i = 0; i < n; i++) {
		w[i] = interval(w1[i].inf(), w2[i].sup());
	}
	return true;
}

bool
//%----------------------------------------------------------------------------
//% Kolev's direct method
//%----------------------------------------------------------------------------
kolevDirect(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w,
	dvector& t1,
	dvector& t3,
	dmatrix& t2,
	dvector& x0)
{
	int n = b1[0].size();
	int K = p1.size();
	omvector M(K + 1);
	ovector b(K + 1), M2(K);
	dvector q(K);
	ivector p(K);
	imatrix B(n, n);
	dmatrix R(n, n), D(n, n), I(n, n), B0(n, K), Hc(n, n), Hr(n, n), B0m(n, K), C(n, n), Y(n, n);

	ktransform(M1, b1, p1, M, b, p);

	R = M[0];
	if (inv(R)) {
		x0 = R * b[0];
		D.fill_in(0.0);
		for (int k = 1; k < K + 1; k++) {
			Y = R * M[k];
			mag(Y, C);
			D = D + C;
			M2[k - 1] = R * (b[k] - M[k] * x0); //% residual correction
		}
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < K; k++) {
				B0(i, k) = M2[k][i];
			}
		}
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) + D(i, j) * interval(-1.0, 1.0); //% B is OK
			}
		}
		if (invRohn(B)) {
			mid(B, Hc);
			rad(B, Hr);
			mag(B0, B0m);
			q.fill_in(1.0);

			t1 = x0;
			t2 = Hc * B0;
			t3 = Hr * B0m * q;
			reduce(t1, t2, t3, w);
			return true;
		}
		else {
			return false;
		}
	}
	return false;
}

void
//%----------------------------------------------------------------------------
//% Another method for computing p-solution it is very sensitive to the 
//% unertainty it seems that it works only for very small unertainty
//%----------------------------------------------------------------------------
kolevL(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w,
	dvector& t1, dvector& t3, dmatrix& t2, dvector& x0)
{
	int n = b1[0].size();
	int niter = 0;
	int K = p1.size();
	real eps = 1.0e-6;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n), I(n, n);

	ivector v(n), v1(n);
	dvector t3i(n);
	dmatrix C0(n, K), C(n, K), L1(n, K), Ap(n, K);
	ovector B(K), B1(K);
	ovector d(K), G(K);
	omvector Q(n), T2(K);

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	ktransform(M1, b1, p1, M, b, p);

	//spectral(M);

	A0 = M[0]; // midpoint matrix of the new system
	b0 = b[0]; // midpoint vector of the new system

	x0 = dvector(n);

	Unit(I);
	R = I - A0; // A(0)

	// v^(0)=x0 - the midpoint solution
	gausse(A0, b0, x0);
	dvector c0 = x0;

	// v^(1) = c^1 + L^1*p
	dvector c1 = R * c0 + b[0];
	for (int k = 0; k < K; k++) {
		B[k] = M[k + 1] * c0 * (-1.0); // B^k=-A^k*c0
	}
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < K; k++) {
			L1(i, k) = B[k][i] + b[k + 1][i];
		}
	}
	v = c1 + L1 * p;

	// v^(2)=t1+t2*p+t3*[-1,1]    
	t1 = dvector(n);
	t3 = dvector(n);
	t2 = dmatrix(n, K);

	t1.fill_in(0.0);
	t2.fill_in(0.0);
	t3.fill_in(0.0);

	t1 = R * c1 + b[0];

	for (int k = 0; k < K; k++) {
		B[k] = M[k + 1] * c1 * (-1.0); // B1^k=-A^k*c1
	}
	dmatrix AL = R * L1;

	for (int i = 0; i < n; i++) {
		for (int k = 0; k < K; k++) {
			t2(i, k) = B[k][i] + b[k + 1][i] + AL(i, k);
		}
	}

	// Allocation of memory for matrix Q
	for (int i = 0; i < n; i++) {
		Q[i] = dmatrix(K, K);
	}

	// Definition of matrix Q
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < K; k++) {
			for (int kk = 0; kk < K; kk++) {
				real s = 0.0;
				for (int j = 0; j < n; j++) {
					s = s + M[k + 1](i, j) * L1(j, kk);
				}
				Q[i](k, kk) = -s;
			}
		}
	}

	// Definition of vector t3
	for (int i = 0; i < n; i++) {
		real s = 0.0;
		for (int j = 0; j < K; j++) {
			for (int k = 0; k < j; k++) {
				s = s + ISLIB_ABS(Q[i](j, k) + Q[i](k, j));
			}
			t3[i] = s;
		}
	}

	// Linearisation
	for (int i = 0; i < n; i++) {
		real s = 0.0;
		for (int k = 0; k < K; k++) {
			t1[i] = t1[i] - 0.5 * Q[i](k, k);
			s = s + 0.5 * fabs(Q[i](k, k));
		}
		t3[i] = t3[i] + s;
	}

	dvector z(n);
	dmatrix R2(n, n);
	dmatrix M2(n, n);

	//============================== Iteration =================================
	reduce(t1, t2, t3, v);

	do {
		niter++;
		v1 = v;

		// Definition of a general B
		for (int k = 0; k < K; k++) {
			B[k] = M[k + 1] * t1 * (-1.0);
		}
		// Calcule the Q values
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < K; k++) {
				for (int kk = 0; kk < K; kk++) {
					real s = 0.0;
					for (int j = 0; j < n; j++) {
						s = s + M[k + 1](i, j) * t2(j, kk);  // Here we use t2 because L2 has other values
					}
					Q[i](k, kk) = -s;
				}
			}
		}

		// Calculation of next C
		t1 = R * t1 + b[0];

		AL = R * t2;

		// Definition of the general L
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < K; k++) {
				t2(i, k) = B[k][i] + AL(i, k) + b[k + 1][i];
			}
		}

		z = t3;
		mag(R, R2);
		t3 = R2 * t3;
		// t3 adds the value of the sum of A*s
		for (int i = 0; i < K; i++) {
			mag(M[i], M2);
			t3 = t3 + M2 * z;
		}

		// t3 adds the values of the linearisation
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < K; j++) {
				for (int k = 0; k < j; k++) {
					t3[i] = t3[i] + ISLIB_ABS(Q[i](j, k) + Q[i](k, j));
				}
			}
		}

		// t3 adds the values of delta (s) and c2 adds the values of linear parts    
		for (int i = 0; i < n; i++) {
			real s = 0.0;
			for (int j = 0; j < K; j++) {
				t1[i] = t1[i] - 0.5 * Q[i](j, j);
				s = s + 0.5 * fabs(Q[i](j, j));
			}
			t3[i] = t3[i] + s;
		}

		reduce(t1, t2, t3, v);

	} while ((stopCrit(v, v1) > eps) && (niter < 5));
	std::cout << "Niter: " << niter << std::endl;

	w = v;

}

void
//%----------------------------------------------------------------------------
//% Original Kolev's Quadratic method. Calculates only outer solution.
//% The system is transformed do the system with parameters varying within 
//% the interval [-1,1]
//%----------------------------------------------------------------------------
kolevQ(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w)
{
	int n = b1[0].size();
	int niter = 0;
	int K = p1.size();
	real eps = 1.0e-6;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n);
	dvector t1, t3, x0;
	dmatrix t2, t4;

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	t1 = dvector(n);
	t3 = dvector(n);
	t2 = dmatrix(n, K);
	t4 = dmatrix(n, K);

	ktransform(M1, b1, p1, M, b, p);
	A0 = M[0]; //% midpoint matrix of the new system
	b0 = b[0]; //% midpoint vector of the new system
	R = A0;
	x0 = dvector(n);
	if (inv(R)) { //% if the mipoint matrix is nonsingular, the computation continues

		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), H(K);
		ovector d(K), G(K);

		gausse(A0, b0, x0); //% computing x0 - the midpoint solution

		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); //% d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; //% B^k=R*A^k
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		//% v^(0)=0
		//% v^(1)
		v = C0 * p;

		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);
		t4.fill_in(0.0);

		//% v^(2)=t1+t2*p+t3
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0; //% wyglada na to, ze D tez jest dobrze obliczone
		}

		//% Calculate t4 (Q) matrix
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < K; k++) {
				t4(i, k) = -D[k](i, k);
			}
		}

		//% t2
		t2 = C0;

		//% Calculate t3 vector (interval part)
		for (int i = 0; i < n; i++) { // n is the size of v
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					s += fabs(D[k](i, j) + D[j](i, k));
				}
			}
			t3[i] = s;
		}
		kolevReduceQ(t1, t2, t3, t4, v);

		do {
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				D[k] = B[k] * t2;
				H[k] = B[k] * t4;
			}

			//% t4 and t2
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < K; k++) {
					t4(i, k) = -D[k](i, k); //% a2- quadratic term
					t2(i, k) = -H[k](i, k) + C0(i, k); //% a3 + a1 - linear term obtained from the linearisation of a3*x^3 + a1*x
				}
			}

			//% t3
			dvector t3i(n);
			for (int i = 0; i < n; i++) {
				real s = 0.0;
				real ss = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < K; j++) {
						if (j != k) {
							s += ISLIB_ABS(H[j](i, k));
						}
					}
				}
				for (int k = 0; k < K; k++) {
					for (int j = k + 1; j < K; j++) {
						s += ISLIB_ABS(D[k](i, j) + D[j](i, k));
					}
				}
				for (int k = 0; k < n; k++) {
					for (int j = 0; j < K; j++) {
						s += ISLIB_ABS(B[j](i, k) * t3[k]); //% tu jest zle
					}
				}
				for (int k = 0; k < K; k++) {
					ss += ISLIB_ABS(H[k](i, k));
				}
				t3i[i] = s + (2.0 * sqrt(3.0) / 9.0) * ss;
			}
			t3 = t3i;

			////////////////////////////////////////////////////////////////////////
			//% ITERATION: computing v^(k), k>=3
			kolevReduceQ(t1, t2, t3, t4, v);
		} while (stopCrit(v, v1) > eps);

		kolevReduceQ(t1, t2, t3, t4, v);
		for (int i = 0; i < n; i++) { //% final result for outer solution
			w[i] = v[i] + interval(x0[i]);
		}
	}
	//std::cout << "Iterations: " << niter << std::endl;
}

void
//%----------------------------------------------------------------------------
//% Original Kolev's Quadratic method
//% Calculates only outer solution
//%----------------------------------------------------------------------------
kolevQ(const omvector& M1, 
	const ovector& b1, 
	const ivector& p1, ivector& w, 
	dvector& t1,
	dvector& t3, 
	dmatrix& t2, 
	dmatrix& t4, 
	dvector& x0)
{
	int n = b1[0].size();
	int niter = 0;
	int K = p1.size();
	real eps = 1.0e-6;
	dvector b0(n);
	dmatrix A0(n, n), R(n, n);

	omvector M = M1;
	ovector b = b1;
	ivector p = p1;

	t1 = dvector(n);
	t3 = dvector(n);
	t2 = dmatrix(n, K);
	t4 = dmatrix(n, K);

	ktransform(M1, b1, p1, M, b, p);
	A0 = M[0]; //% midpoint matrix of the new system
	b0 = b[0]; //% midpoint vector of the new system
	R = A0;
	x0 = dvector(n);
	if (inv(R)) { //% if the mipoint matrix is nonsingular, the computation continues

		ivector v(n), v1(n);
		dvector t3i(n);
		dmatrix C0(n, K), C(n, K);
		omvector B(K), D(K), H(K);
		ovector d(K), G(K);

		gausse(A0, b0, x0); //% computing x0 - the midpoint solution

		for (int k = 0; k < K; k++) {
			d[k] = R * (b[k + 1] - M[k + 1] * x0); //% d^k = R*(b^k-A^k*x0)
			B[k] = R * M[k + 1]; //% B^k=R*A^k
			for (int j = 0; j < n; j++) {
				C0(j, k) = d[k][j];
			}
		}
		//% v^(0)=0
		//% v^(1)
		v = C0 * p;

		t1.fill_in(0.0);
		t2.fill_in(0.0);
		t3.fill_in(0.0);
		t4.fill_in(0.0);
		//% v^(2)=t1+t2*p+t3
		for (int k = 0; k < K; k++) {
			D[k] = B[k] * C0; // wyglada na to, ze D tez jest dobrze obliczone
		}
		//% Calculate t4 (Q) matrix
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < K; k++) {
				t4(i, k) = -D[k](i, k);
			}
		}
		//% t2
		t2 = C0;
		//% Calculate t3 vector (interval part)
		for (int i = 0; i < n; i++) { //% n is the size of v
			real s = 0.0;
			for (int k = 0; k < K; k++) {
				for (int j = k + 1; j < K; j++) {
					s += fabs(D[k](i, j) + D[j](i, k));
				}
			}
			t3[i] = s;
		}
		kolevReduceQ(t1, t2, t3, t4, v);

		do {
			niter++;
			v1 = v;
			for (int k = 0; k < K; k++) {
				D[k] = B[k] * t2;
				H[k] = B[k] * t4;
			}

			//% t4 and t2
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < K; k++) {
					t4(i, k) = -D[k](i, k); //% a2- quadratic term
					t2(i, k) = -H[k](i, k) + C0(i, k); //% a3 + a1 - linear term obtained from the linearisation of a3*x^3 + a1*x
				}
			}

			//% t3
			dvector t3i(n);
			for (int i = 0; i < n; i++) {
				real s = 0.0;
				real ss = 0.0;
				for (int k = 0; k < K; k++) {
					for (int j = 0; j < K; j++) {
						if (j != k) {
							s += ISLIB_ABS(H[j](i, k));
						}
					}
				}
				for (int k = 0; k < K; k++) {
					for (int j = k + 1; j < K; j++) {
						s += ISLIB_ABS(D[k](i, j) + D[j](i, k));
					}
				}
				for (int k = 0; k < n; k++) {
					for (int j = 0; j < K; j++) {
						s += ISLIB_ABS(B[j](i, k) * t3[k]); //% tu jest zle
					}
				}
				for (int k = 0; k < K; k++) {
					ss += ISLIB_ABS(H[k](i, k));
				}
				t3i[i] = s + (2.0 * sqrt(3.0) / 9.0) * ss;
			}
			t3 = t3i;

			////////////////////////////////////////////////////////////////////////
			// ITERATION: computing v^(k), k>=3
			kolevReduceQ(t1, t2, t3, t4, v);
		} while (stopCrit(v, v1) > eps);

		kolevReduceQ(t1, t2, t3, t4, v);
		for (int i = 0; i < n; i++) { //% final result for outer solution
			w[i] = v[i] + interval(x0[i]);
		}
	}
	//std::cout << "Iterations: " << niter << std::endl;
} //% >>> kolevQ <<<