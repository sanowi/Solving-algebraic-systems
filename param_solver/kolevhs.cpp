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
#include "../evolution/evolution.h"
#include "../param_solver/monotonicity.h"
#include "kolev.h"
#include "kolevhs.h"

#define EPSP 1.0e-6
#define EPSQ 1.0e-6
#define MAXITER 1000

#define EVOL 0
#define tau 0.05

void eIntHS(int minmax, int m, const dvector& x0, const dvector& t1, const dmatrix& t2, const dvector& t3, interval& eI, interval& mw)
//%
//% Computes the interval which encloses left (minmax==-1)/right (minmax==1) 
//% endpoint of the hull
//%
{
	int n = x0.size();
	int K = t2.num_cols();
	real t2m = 0.0;

	eI = x0[m] + t1[m];
	for (int i = 0; i < K; i++) {
		t2m = t2m + ISLIB_ABS(t2(m, i)); // +|Lj|
	}
	mw = eI + (t2m + t3[m]) * interval(-1.0, 1.0);
	eI = eI + minmax * t2m + t3[m] * interval(-1.0, 1.0);
}

real chkSolHS(int minmax, int m, const interval& eI, const dmatrix& t2, const ivector& pi, const omvector& M1, const ovector& b1)
{
	int n = t2.num_rows();
	int K = pi.size();
	dvector bs, psol(K), ws(n);
	dmatrix As;
	if (minmax == -1) {
		for (int k = 0; k < K; k++) {
			if (t2(m,k) > 0) {
				psol[k] = pi[k].inf();
			}
			else {
				psol[k] = pi[k].sup();
			}
		}
	}
	else {
		for (int k = 0; k < K; k++) {
			if (t2(m, k) > 0) {
				psol[k] = pi[k].sup();
			}
			else {
				psol[k] = pi[k].inf();
			}
		}
	}
	omvec2real(M1, psol, As);
	ovec2real(b1, psol, bs);
	gausse(As, bs, ws);
	return ws[m];
}

void
contractorHS(real& red, int minmax, int m, ivector& p1, const dvector& innl, const dvector& innu, const dvector& t1, const dvector& t3,
const dmatrix& t2, const dvector& x0)
//%
//% Kolev's version slightly modified
//% the m-th entry of the solution is optimized
//%
{
	int i, n = t1.size(), K = p1.size();
	real t2sum = 0.0;
	interval d;
	dvector ha(K);
	ivector p(K);
	p.fill_in(interval(-1.0, 1.0));

#if (EVOL == 0) 
	//%
	//% computing bounds on the hull using
	//% the the lower bound of the inner solution and the lower bound of the outer solution for left endpoint
	//% and the the upper bound of the inner solution and the upper bound of the outer solution for right endpoint
	//% 
	interval v;
	kolevReducei(m, t1, t2, t3, v);
	if (minmax == -1) {
		//% upper bound for the hull solution <- left endpoint of thr inner solution
		d = interval(v.inf() - t1[m], innl[m] - t1[m] - x0[m]); 
	}
	else {
		//% lower bound of the hull solution <- right endpoint of the inner solution
		d = interval(innu[m] - t1[m] - x0[m], v.sup() - t1[m]); 
	}
#elif (EVOL == 1) 
	//%
	//% computing bound on hull using evolutionary algorithm
	//%
	real inf, d0;
	omatrix oA(n, n);
	ovector ob(n);
	ivector q(1); q[0] = 0.0;
	interval v;
	kolevReducei(m, t1, t2, t3, v);
	for (i = m; i < m + 1; i++) {
		evolution(minmax, i, p1, q, oA, ob, inf);
	}
	d0 = minmax * inf - (x0[m] + t1[m]);
	if (minmax == -1) {
		d = interval(v.inf() - t1[m], d0);
	}
	else {
		d = interval(d0, v.sup() - t1[m]);
	}
#else 
	//%
	//%  Kolev's version: computing upper bound 
	//% using endpoint of outer solution modified by tau*t3[m]
	//%
	interval v;
	kolevReducei(m, t1, t2, t3, v);
	if (minmax == -1) {
		d = interval(v.inf() - t1[m] + tau * t3[m]); // upper bound for the hull solution <- left endpoint of thr inner solution
	}
	else {
		d = interval(v.sup() - t1[m] - tau * t3[m]); // lower bound for the hull solution <- right endpoint of the inner solution
	}
#endif
	for (i = 0; i < K; i++) ha[i] = ISLIB_ABS(t2(m, i)); // at2m = |L_m:| = |t2(m,:)| <-> ha
	int kmax = 0;
	for (i = 1; i < K; i++) if (ha[i] > ha[kmax]) kmax = i;

	interval sum_t2p = 0.0;
	for (i = 0; i < K; i++) {
		if (i != kmax) {
			sum_t2p = sum_t2p + t2(m, i) * p[i]; // L<->t2
		}
	}
	//%-------------------------------------------------------------------------------
	interval bb = d + sum_t2p + interval(-t3[m], t3[m]); // in Kolev [b]
	interval qq = bb / t2(m, kmax);
	//%-------------------------------------------------------------------------------

	//% After reduction of u=[-1,1], the original parameter must be reduced accordingly
	interval pk = interval(ISLIB_MAX(qq.inf(), p[kmax].inf()), ISLIB_MIN(qq.sup(), p[kmax].sup()));
	real pkd = pk.inf();
	real pkg = pk.sup();
	real rk = pk.rad(); // rk
	//% Case 1
	if ((pkd>-1.0) && (pkg < 1.0)) { // both sides reduction, pk strictly included in p1
		std::cout << "Both sides reduction" << std::endl;
		red = 1.0 - rk;
		p1[kmax] = interval(p1[kmax].inf() + p1[kmax].rad()*(pkd + 1.0), p1[kmax].inf() + p1[kmax].rad()*(pkg + 1.0));;
		return;
	}
	//% Case 2
	else if ((pkd > 1.0) || (pkg < -1.0)) {
		// no intersection
		std::cout << "No intersection" << std::endl;
		red = -200;
		return;
	}
	//% Case 3
	else {
		real rpmax = p1[kmax].rad(); // r_k
		real rs = rpmax * rk;
		real pc;
		red = 1.0 - rk;
		if (pkd > -1.0) { // left reduction
			//std::cout << "left side reduction: ";
			pc = p1[kmax].sup() - rs; 
		}
		else if (pkg < 1.0) { // right reduction
			//std::cout << "right side reduction: ";
			pc = p1[kmax].inf() + rs; 
		}
		else {
			// No reduction
			red = 0;
			return;
		}
		p1[kmax] = interval(pc - rs, pc + rs);
	}
}

void 
kolevHS(int m, int minmax, const omvector& M1, const ovector& b1, const ivector& p1, interval& mw, dvector& pd, interval& eI, int& iter)
// calculates m-th entry of the hull solution
{
	int n = b1[0].size(), K = p1.size(), niter;
	real xdm, xgm, maxr, red;
	dvector b0(n), t1(n), t3(n), x0(n);
	dvector pp(K), bi(n), xa(n), innl(n), innu(n);
	dmatrix A0(n, n), A(n, n), t2(n, K);

	ivector p(K), pi(K), hx(n), w(n), iw(n), ww(n);
	omvector M(K + 1), Mi(K + 1);
	ovector b(K + 1);

	pi = p1;
	for (iter = 0; iter < MAXITER; iter++) { // loops until q([p],[p]')<=eps_q
		p = pi;
		//--------------------------------------------------------------------------
		bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter); // outer p-solution
		//--------------------------------------------------------------------------
		kolevInner(x0, t1, t3, t2, innl, innu); // innl, innu - left and right endpoint of the inner solution
		//--------------------------------------------------------------------------
		xdm = w[m].inf(); xgm = w[m].sup(); // xdm corr. to xdi0, xgm corr. xgi0
		maxr = maxrad(pi);
		if (maxr < EPSP) {
			std::cout << "SUCCESS!!! in " << iter << " iterations" << std::endl;
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			std::cout << "pi" << std::endl << pi << std::endl;
			eIntHS(minmax, m, x0, t1, t2, t3, eI, mw);
			real mw = chkSolHS(minmax, m, eI, t2, pi, M1, b1);
			//std::cout << "eI = " << eI << std::endl;
			//std::cout << "mw = " << mw << std::endl;
			if (mw >= eI.inf() && mw <= eI.sup()) {
				std::cout << "RIGHT!!!" << std::endl;
			}
			else {
				std::cout << "WRONG!!!: error = " << ISLIB_ABS(mw - eI.mid()) << std::endl;
			}
			return;
		}
		contractorHS(red, minmax, m, pi, innl, innu, t1, t3, t2, x0);
		//std::cout << "red = " << red << std::endl;
		if (red < EPSQ) { //% no improvement
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << "CRUDE!!!: in " << iter << " iterations" << std::endl;
			std::cout << "pi" << std::endl << pi << std::endl;
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			mw = w[m];
			std::cout << "mw: " << mw << std::endl;
			Monot(pi, M1, b1, ww);
			std::cout << "ww" << std::endl << ww << std::endl;
			std::cout << "ww[m]: " << ww[m] << std::endl;
			mw = interval(ISLIB_MAX(mw.inf(), ww[m].inf()), ISLIB_MIN(mw.sup(),ww[m].sup()));
			std::cout << "mw: " << mw << std::endl;
			// Tu mozna jeszcze probowc liczyc tak jak to Kolev robil
			// Trzeba napisac kodzik
			std::cout << "----------------------------------------------" << std::endl;
			return;
		}
		//char c = _getch();
	}
	std::cout << "FAILURE: maximum number of iterations exceeded!!!" << std::endl;
}

//////////////////////////////////////////////////
// MOJA WERSJA - NA RAZIE NIE DZIALA POPRAWNIE
//////////////////////////////////////////////////
void
//% contracts the parameters, corresponds to the part of the procedure P2
//% the matrix oA and the vector ob must be formes automatically from M1 and b1!!!!
//% m - number of the entry to be optmized
//% p1 - vector of original interval parameters, it is REQUIRED(!!!) for evolutionary algorithm
//% p - normalized interval parameters (denoted by u in Kolev's paper (ACM))
//% constraint equation for m-th entry of the parametric solution has the form
//% \sum_{j=1}^m l_{kj}p_j + s_k = d_k, where s_k is the same as t3[k] vector in this program
contractor(int minmax, 
	int m, 
	ivector& p, 
	const ivector& p1, 
	const ivector& ww, 
	const dvector& t1, 
	const dvector& t3,
	const dmatrix& t2, 
	const dvector& x0, 
	int fun)
{
	int i, n = ww.size(), K = p.size(); // we optimize the m-th entry of the solution
	interval ptmp;
	ivector q(1); q[0] = 0.0;
	real inf;
	omatrix oA(n, n);
	ovector ob(n);

	// computing upper (minmax=-1) or lower (minmax=1) bound on the m-th entry of the hull
	for (i = m; i < m + 1; i++) {
		evolution(minmax, i, p1, q, oA, ob, inf, fun);
	}
	real d = minmax * inf - (x0[m] + t1[m]); // d_k in Kolev's paper with evolution
	int k; // k indicates h with the greatest absolute value
	for (k = 0; k < K; k++) {
		interval sum_hip = 0.0;
		for (i = 0; i < K; i++) {
			if (i != k) {
				sum_hip = sum_hip + t2(m, i) * p[i]; // t2 corresponds to L in the paper
			}
		}
		interval bb = d - sum_hip - interval(-t3[m], t3[m]); // u Koleva po prostu [b]
		interval qq = bb / t2(m, k);
		ptmp = interval(ISLIB_MAX(qq.inf(), p[k].inf()), ISLIB_MIN(qq.sup(), p[k].sup()));
		// for kolevSWIM method, this must be commented
		if (!ptmp.is_empty()) {
			p[k] = ptmp;
		}
	}
}

void kolevSWIM(int minmax, 
	const omvector& M1, 
	const ovector& b1, 
	const ivector& p1, 
	ivector& w,
	int fun)
// calculates i-th entry of the hull solution
{
	int n = b1[0].size(), K = p1.size(), iter = 0, niter;
	dvector t1(n), t3(n);
	dmatrix t2(n, K);

	dvector b0(n), x0(n);
	dmatrix A0(n, n);
	ivector p(K), pi(K), hx(n), iw(n);
	omvector M(K + 1), Mi(K + 1);
	ovector b(K + 1), bi(K + 1);
	interval eiI;

	omatrix oA(n, n);
	ovector ob(n);

	w = ivector(n);
	hx.fill_in(0.0);
	for (int m = 0; m < n; m++) { // loops over entries of the solution
		//std::cout << "===========================" << std::endl << "m = " << m << std::endl;
		//std::cout << "===========================" << std::endl;
		// since each entry is optimized separately, we have to start
		// computation for each entry with initial  values of the system matrix,
		// right-hand vector and parameter vector
		pi = p1;
		Mi = M1;
		bi = b1;
		//do {
		//	iter++;
		for (iter = 0; iter < MAXITER; iter++) { // loops until q([p],[p]')<=eps_q
			// FIRST STEP: the system is transformed so that all [p]=[-1,1]
			ktransform(Mi, bi, pi, M, b, p);
			// SECOND STEP: the system is solved using Kolev's method
			// we obtain a parametric solution of the form
			// [x](p)=[l](p)=x0+t1+t2*p+[-t3,t3], p\in[p]
			bool conv = kolevOuter(M, b, p, w, t1, t3, t2, x0, niter);
			//std::cout << "w(kolevS)" << std::endl << w << std::endl;
			// this can be removed late since we use the result of EOM as a lower (upper) bound
			// in the constraint equation (see THIRD STEP)
			ivector ww(n);
			ww.fill_in(0.0);
			kolevInner(x0, t1, t3, t2, iw);
			if (minmax = -1) {
				eiI = interval(w[m].inf(), iw[m].inf());
			}
			else {
				eiI = interval(iw[m].sup(), iw[m].sup());
			}
			std::cout << "eiI = " << eiI << std::endl;
			//kolevInner(x0, t1, t3, t2, ww);
			// THIRD STEP: contract parameter intervals using 
			// the constraint equation for k-th entry of the p-solution
			// \sum_{j=1}^m L_{kj}+[\rho_k,\rho_k]=d_k,
			// we have to do it separately for each entry of the solution and for upper and lower bound
			// where d_k=l^u_k-(x^(0)_k+c_k), and l^u_k is obtained from the evolutionary algorithm
			// first parameter of the contracting method indicates whether 
			// lower (-1) or upper (+1) 
			// bound on the hull entry is computed
			// p1 - vector of original parameters is REQUIRED for EOM
			// p - vector of intervals [-1,1]
			contractor(minmax, m, p, p1, ww, t1, t3, t2, x0, fun);
			hx[m] = w[m];
			//std::cout << "(!!!)p = " << std::endl << p << std::endl << std::endl;
			//std::cout << "maxrad = " << maxrad(p) << std::endl;
			if (maxrad(p) < EPSP) {
				// bound for m-th entry found
				hx[m] = w[m];
				//if (minmax = -1) hx[m] = interval(w[m].inf());
				//if (minmax = 1) hx[m] = interval(w[m].sup());
				std::cout << "hx[m] = " << w[m] << std::endl;
				break;
			}
			pi = p;
			Mi = M;
			bi = b;
		} //while (maxrad(p) > EPSP & iter < 200);
	}
	std::cout << "Aaaaa" << std::endl;
	//std::cout << "hx[m] = " << w[m] << std::endl;
	w = hx; //std::cout << hx << std::endl;
}

void testContractor(int m, const dmatrix& t2, const dvector& t3, const ivector& p, const interval& d)
{
	int i, n = t2.num_rows(), K = t2.num_cols();
	dvector ha(K);

	for (i = 0; i < K; i++) ha[i] = ISLIB_ABS(t2(m, i)); // at2m = |L_m:| = |t2(m,:)| <-> ha
	int kmax = 0;
	for (i = 1; i < K; i++) if (ha[i] > ha[kmax]) kmax = i;
	std::cout << "kmax = " << kmax << std::endl;

	interval sum_t2p = 0.0;
	for (i = 0; i < K; i++) {
		for (int j = 0; j < K; j++) {
			if (i != j) {
				sum_t2p = sum_t2p + t2(m, i) * p[i]; // L<->t2
			}
		}
		interval bb = d + sum_t2p + interval(-t3[m], t3[m]); // in Kolev [b]
		std::cout << "bk (" << i << ") = " << bb << ", Rk = " << (sum_t2p + interval(-t3[m], t3[m])).mag() << std::endl;
		interval qq = bb / t2(m, i);
		std::cout << "bf(" << i << ") after = " << qq << std::endl;
	}
	char c = _getch();
	//%-------------------------------------------------------------------------------
}