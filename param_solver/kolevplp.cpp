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
#include "../param_solver/iterative_psolvers.h"
#include "kolev.h"
#include "kolevplp.h"

#define EPSP 1.0e-6
#define EPSQ 1.0e-6
#define MAXITER 1000

#define EVOLPLP 0
#define tau1 0.0 //right // for rho=0.1,...,0.263 tau1=0.18, tau2=0.02
#define tau2 0.0 //left
#define tau11 0.5 //4 //0.9
#define tau22 0.5 //02 // 0.02

//% For example ex3 (exSWIM) uncertainty is: max PLP = 0.269, max PLP(Kolev bound) = 0.293

void
eIntPLP(int minmax, const dvector& c, const dvector& x0, const dvector& t1, const dmatrix& t2, const dvector& t3, interval& eI)
//%
//% Computes the interval which encloses left (minmax==-1)/right (minmax==1) 
//% endpoint of the hull
//%
{
	int n = c.size();
	int K = t2.num_cols();
	real t3a = 0.0;
	real at2ma = 0.0;
	dvector t2m(K);

	for (int k = 0; k < K; k++) {
		real s = 0.0;
		for (int i = 0; i < n; i++) {
			s = s + c[i] * t2(i, k);
		}
		t2m[k] = s;
		at2ma = at2ma + ISLIB_ABS(t2m[k]); // +|Lj|
	}
	eI = 0.0;
	for (int i = 0; i < n; i++) {
		eI = eI + c[i] * (x0[i] + t1[i]);
		t3a = t3a + ISLIB_ABS(c[i]) * t3[i];
	}
	eI = eI + (at2ma + t3a) * interval(-1.0, 1.0); // interval that encloses lower/upper endpoint of the hull
}

void
contractorPLP_Working(real& red, int minmax, const dvector& c, ivector& p1, const dvector& innl, 
const dvector& innu, const dvector& t1, const dvector& t3, const dmatrix& t2, const dvector& x0)
//%
//% Kolev's version slightly modified
//%
{
	int i, n = t1.size(), K = p1.size(); // we optimize the m-th entry of the solution
	real t2sum = 0.0;
	interval d;
	dvector at2m(K), t2m(K);
	ivector p(K);

	omatrix oA(n, n);
	ovector ob(n);
	ivector q(1); q[0] = 0.0;

	p.fill_in(interval(-1.0, 1.0));

	for (i = 0; i < K; i++) {
		t2m[i] = 0.0;
		for (int m = 0; m < n; m++) {
			t2m[i] += c[m] * t2(m, i); // Lj
		}
		at2m[i] = ISLIB_ABS(t2m[i]); // |Lj|
	}
	real t3m = 0.0;
	for (int m = 0; m < n; m++){
		t3m = t3m + ISLIB_ABS(c[m]) * t3[m]; // 
	}
	interval xi0;
	real d0 = 0.0;
	for (int i = 0; i < n; i++) {
		d0 = d0 + c[i] * (x0[i] + t1[i]);
	}
	xi0 = d0;
	for (int i = 0; i < K; i++) {
		xi0 = xi0 + at2m[i] * interval(-1.0, 1.0);
	}
	for (int i = 0; i < n; i++) {
		xi0 = xi0 + t3m * interval(-1.0,1.0);
	}

#if (EVOLPLP == 1)
	real inf;
	evolution(minmax, 0, p1, q, oA, ob, inf);
	if (minmax==-1) {
		d = interval(xi0.inf() - d0, minmax * inf - d0);
	}
	else {
		d = interval(minmax * inf - d0, xi0.sup() - d0);
	}
#elif (EVOLPLP == 2)
	if (minmax == -1) {
		d = t3m * interval(-1.0, 1.0);
		for (int m = 0; m < K; m++) {
			d = d - at2m[m];
		}
	}
	else {
		d = t3m * interval(-1.0, 1.0);
		for (int m = 0; m < K; m++) {
			d = d + at2m[m];
		}
	}
#elif (EVOLPLP == 3)
	real inf;
	evolution(minmax, 0, p1, q, oA, ob, inf);
	//d = minmax * inf - d0; //(x0[m] + t1[m]);// tu trzeba odejmowac somu
	if (minmax==-1) {
		//d = interval(xi0.inf() - d0, minmax * inf - d0 + 0.002 * t3m); // 0.002 for 0.32 ex3
		//d = interval(xi0.inf() - d0, minmax * inf - d0 - 0.05 * t3m); // 0.33 ex3
		d = interval(xi0.inf() - d0, minmax * inf - d0 - tau11 * t3m); // 0.34 ex3
	}
	else {
		//d = interval(minmax * inf - d0 + tau11 * t3m, xi0.sup() - d0);
		d = interval(minmax * inf - d0 + tau22 * t3m, xi0.sup() - d0); // 0.32, 0.33 ex3
		//d = interval(minmax * inf - d0 - 0.05 * t3m, xi0.sup() - d0);
	}
#else // Kolev's constraint propagation version
	if (minmax == -1) {
		d = tau1 * t3m;
		for (int m = 0; m < K; m++) {
			d = d - at2m[m];
		}
	}
	else {
		d = -tau2 * t3m;
		for (int m = 0; m < K; m++) {
			d = d + at2m[m];
		}
	}
#endif
	int kmax = 0;
	for (i = 1; i < K; i++) if (at2m[i] > at2m[kmax]) kmax = i;

	interval sum_t2p = 0.0;
	for (i = 0; i < K; i++) {
		if (i != kmax) {
			sum_t2p = sum_t2p + t2m[i] * p[i]; // c*L*p, t2 corresponds to L in the paper
		}
	}
	//-------------------------------------------------------------------------------

	interval bb = d - sum_t2p - interval(-t3m, t3m); // u Koleva po prostu [b]
	//std::cout << "bb = " << bb << std::endl;
	interval qq = bb / t2m[kmax];
	//std::cout << "qq = " << qq << std::endl;
	//-------------------------------------------------------------------------------

	interval pk = interval(ISLIB_MAX(qq.inf(), p[kmax].inf()), ISLIB_MIN(qq.sup(), p[kmax].sup()));
	real bkd = bb.inf();
	real bkg = bb.sup();
	real pkd = pk.inf();
	real pkg = pk.sup();
	real rpk = pk.rad(); // rk

	if ((pkd>-1.0) && (pkg < 1.0)) {
		// pKI strictly included in p1
		red = 1 - rpk;
		p1[kmax] = interval(p1[kmax].inf() + p1[kmax].rad()*(pkd + 1.0), p1[kmax].inf() + p1[kmax].rad()*(pkg + 1.0));;
		return;
	}
	else if ((pkd > 1.0) || (pkg < -1.0)) {
		// no intersection
		std::cout << "No intersection!!!" << std::endl;
		red = -200;
		return;
	}
	else {
		real rpmax = p1[kmax].rad(); // r_k
		real rs = rpmax * rpk;
		real pc;
		red = 1.0 - rpk;
		if (pkd > -1.0) { // Left reduction
			pc = p1[kmax].sup() - rs;
		}
		else if (pkg < 1.0) { // Right reduction
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

//%=======================================================================
//% KOLEVOPT is used to decide which bound is used
#define KOLEVOPT 4 //1 - B1, 2 - B2, 3 - B3, 4 - B4(=modf. B3)
//%=======================================================================

void
contractorPLP(real& red, int minmax, const omvector& M, const ovector& b, const dvector& c, ivector& p1, 
const dvector& innl, const dvector& innu, const dvector& t1, const dvector& t3, const dmatrix& t2, const dvector& x0)
//%
//% Kolev's version slightly modified
//%
{
	int i, n = t1.size(), K = p1.size(); // we optimize the m-th entry of the solution
	real t2sum = 0.0;
	interval d; // d is like l^u in the paper
	dvector at2m(K), t2m(K);
	ivector p(K);

	omatrix oA(n, n);
	ovector ob(n);
	ivector q(1); q[0] = 0.0;

	p.fill_in(interval(-1.0, 1.0));

	for (i = 0; i < K; i++) {
		t2m[i] = 0.0;
		for (int m = 0; m < n; m++) {
			t2m[i] += c[m] * t2(m, i); // Lj
		}
		at2m[i] = ISLIB_ABS(t2m[i]); // |Lj|
	}
	real t3m = 0.0;
	for (int m = 0; m < n; m++){
		t3m = t3m + ISLIB_ABS(c[m]) * t3[m]; // 
	}
	interval xi0; // xi0 is an outer solution
	real d0 = 0.0;
	for (int i = 0; i < n; i++) {
		d0 = d0 + c[i] * (x0[i] + t1[i]);
	}
	xi0 = d0;
	for (int i = 0; i < K; i++) {
		xi0 = xi0 + at2m[i] * interval(-1.0, 1.0);
	}
	for (int i = 0; i < n; i++) {
		xi0 = xi0 + t3m * interval(-1.0, 1.0);
	}

#if (KOLEVOPT == 1) // Case B1
	if (minmax == -1) {
		d = t3m;
		for (int m = 0; m < K; m++) {
			d = d - at2m[m];
		}
	}
	else {
		d = -t3m;
		for (int m = 0; m < K; m++) {
			d = d + at2m[m];
		}
	}
#elif (KOLEVOPT == 2) // Case B2
	real inf;
	evolution(minmax, 0, p1, q, oA, ob, inf);
	if (minmax == -1) {
		//d = interval(xi0.inf() - d0, minmax * inf - d0);
		d = minmax * inf - d0;
	}
	else {
		//d = interval(minmax * inf - d0, xi0.sup() - d0);
		d = minmax * inf - d0;
	}
#elif (KOLEVOPT == 3)
	ivector pt(K), ww(n);
	if (minmax == 1) { // maximum
		for (int i = 0; i < K; i++) {
			if (t2m[i] > 0) {
				pt[i] = p1[i].sup();
			}
			else {
				pt[i] = p1[i].inf();
			}
		}
	}
	else {
		for (int i = 0; i < K; i++) {
			if (t2m[i] < 0) {
				pt[i] = p1[i].sup();
			}
			else {
				pt[i] = p1[i].inf();
			}
		}
	}
	kolevOuter(M, b, pt, ww);
	d = 0.0;
	for (int i = 0; i < n; i++) {
		d += c[i] * ww[i];
	}
	d = d - d0;
#else // B4

	// Example 3
	//real tauL = 0.74;
	//real tauU = 0.71; // 0.71 is good for 0.3, but not good for 0.4

	real tauL = 0.8;
	real tauU = 0.94;
	ivector pt(K), ww(n);
	if (minmax == 1) { // maximum
		for (int i = 0; i < K; i++) {
			if (t2m[i] > 0) {
				pt[i] = p1[i].sup();
			}
			else {
				pt[i] = p1[i].inf();
			}
		}
	}
	else {
		for (int i = 0; i < K; i++) {
			if (t2m[i] < 0) {
				pt[i] = p1[i].sup();
			}
			else {
				pt[i] = p1[i].inf();
			}
		}
	}
	kolevOuter(M, b, pt, ww);
	d = 0.0;
	for (int i = 0; i < n; i++) {
		d += c[i] * ww[i];
	}
	if (minmax == 1) {
		real ll = xi0.sup();
		d = ll + tauU * (d - ll);
	}
	else {
		real ll = xi0.inf();
		d = ll + tauL * (d - ll);
	}
	d = d - d0;
#endif
	int kmax = 0;
	for (i = 1; i < K; i++) if (at2m[i] > at2m[kmax]) kmax = i;

	interval sum_t2p = 0.0;
	for (i = 0; i < K; i++) {
		if (i != kmax) {
			sum_t2p = sum_t2p + t2m[i] * p[i]; // c*L*p, t2 corresponds to L in the paper
		}
	}
	//-------------------------------------------------------------------------------
	interval bb = d - sum_t2p - interval(-t3m, t3m); // u Koleva po prostu [b]
	//std::cout << t2m << std::endl;
	interval qq = bb / t2m[kmax];
	//-------------------------------------------------------------------------------

	interval pk = interval(ISLIB_MAX(qq.inf(), p[kmax].inf()), ISLIB_MIN(qq.sup(), p[kmax].sup()));
	real bkd = bb.inf();
	real bkg = bb.sup();
	real pkd = pk.inf();
	real pkg = pk.sup();
	real rpk = pk.rad(); // rk

	if ((pkd>-1.0) && (pkg < 1.0)) {
		// pKI strictly included in p1
		red = 1 - rpk;
		p1[kmax] = interval(p1[kmax].inf() + p1[kmax].rad()*(pkd + 1.0), p1[kmax].inf() + p1[kmax].rad()*(pkg + 1.0));;
		return;
	}
	else if ((pkd > 1.0) || (pkg < -1.0)) {
		// no intersection
		std::cout << "No intersection!!!" << std::endl;
		red = -200;
		return;
	}
	else {
		real rpmax = p1[kmax].rad(); // r_k
		real rs = rpmax * rpk;
		real pc;
		red = 1.0 - rpk;
		if (pkd > -1.0) { // Left reduction
			pc = p1[kmax].sup() - rs;
		}
		else if (pkg < 1.0) { // Right reduction
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

#define DIRECT 0

void 
kolevPLP_original_working(int minmax, const dvector& c, const omvector& M1, const ovector& b1, const ivector& p1, 
interval& mw, int& iter, real& time)
//% PLP solver
//% calculates m-th entry of the hull solution
{
	int n = b1[0].size(), K = p1.size(), niter;
	interval eI;
	real maxr, red;

	dvector t1(n), t3(n), pp(K), bi(n), xa(n), innl(n), innu(n), b0(n), x0(n);
	dmatrix t2(n, K), A0(n, n), A(n, n);

	ivector p(K), pi(K), hx(n), w(n), iw(n), ww(n);
	omvector M(K + 1), Mi(K + 1);
	ovector b(K + 1), ob(n);
	omatrix oA(n, n);

	time_t start = clock();

	pi = p1;
	for (iter = 0; iter < MAXITER; iter++) { // loops until q([p],[p]')<=eps_q
		p = pi;
		//--------------------------------------------------------------------------
#if DIRECT == 1
		bool conv = kolevDirect(M1, b1, pi, w, t1, t3, t2, x0);
#else
		bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter); // outer p-solution
#endif
		//--------------------------------------------------------------------------
		// NOTE: we have to take innl and innu since the inner solution can be empty
		kolevInner(x0, t1, t3, t2, innl, innu); // inner solution to get bounds on the hull
		//--------------------------------------------------------------------------
		//xdm = w[m].inf(); xgm = w[m].sup(); // xdm corrsp. to xdi0, xgm corrsp. xgi0
		maxr = maxrad(pi);
		if (maxr < EPSP) {
			std::cout << "SUCCESS!!! in " << iter << " iterations" << std::endl;

			//% computing the final enclosure
#if DIRECT == 1
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
#else
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter); // outer p-solution
#endif

			//% computing the points at which the extrema are attained
			ivector pmid(K);
			for (int kk = 0; kk < K; kk++) {
				//std::cout << pi[kk].inf() << " ";
				std::cout << pi[kk] << " ";
				pmid[kk] = pi[kk].mid();
			}
			//for (int kk = 0; kk < K; kk++) {
			//	//std::cout << pi[kk].sup() << " ";
			//	std::cout << pi[kk] << " ";
			//}
			std::cout << std::endl;
			kolevOuter(M1, b1, pmid, w, t1, t3, t2, x0, niter);
			// to obtain the optimal value range, we can do the following two lines
			/*kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			std::cout << "w" << std::endl << w << std::endl;*/
			if (minmax == 1) {
				std::cout << "max: " << w << std::endl;
			} 
			else {
				std::cout << "min: " << w << std::endl;
			}
			
			//% computing the final solution
			eIntPLP(minmax, c, x0, t1, t2, t3, mw);
			std::cout << "eI = " << mw << std::endl;
			/*if (mws >= eI.inf() && mws <= eI.sup()) {
				std::cout << "RIGHT!!!" << std::endl;
			}
			else {
				std::cout << "WRONG!!!: error = " << ISLIB_ABS(mws - eI.mid()) << std::endl;
			}*/
			time = difftime(clock(), start) / CLOCKS_PER_SEC;
			return;
		}
		//%
		//% Contracting the interval using simple constraint propagation method
		contractorPLP(red, minmax, M1, b1, c, pi, innl, innu, t1, t3, t2, x0);
		//%
		if (red < EPSQ) {
			//%
			//% No improvement: only crude enclosure found
			//% HERE we can try to check the monotonicity 
			std::cout << "CRUDE!!! in " << iter << " iterations" << std::endl;
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].inf() << " ";
			}
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].sup() << " ";
			}
			std::cout << std::endl;
			/*bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			eIntPLP(minmax, c, x0, t1, t2, t3, mw);*/
			
			//std::cout << "pi" << std::endl << pi << std::endl;
			
			//% trying monotonicity
			MonotPLP(minmax, c, pi, M1, b1, mw);
			//std::cout << "CRUDE mmw(monot)" << mw << std::endl;
			return;
		}
	}
	std::cout << "FAILURE: maximum iterations exceeded!!!" << std::endl;
	return;
}

void
kolevPLP(int minmax, const dvector& c, const omvector& M1, const ovector& b1, const ivector& p1,
interval& mw, int& iter, real& time)
//% PLP solver
//% calculates m-th entry of the hull solution
{
	int n = b1[0].size(), K = p1.size(), niter;
	interval eI;
	real maxr, red;

	dvector t1(n), t3(n), pp(K), bi(n), xa(n), innl(n), innu(n), b0(n), x0(n);
	dmatrix t2(n, K), A0(n, n), A(n, n);

	ivector p(K), pi(K), hx(n), w(n), iw(n), ww(n);
	omvector M(K + 1), Mi(K + 1);
	ovector b(K + 1), ob(n);
	omatrix oA(n, n);

	time_t start = clock();

	pi = p1;
	for (iter = 0; iter < MAXITER; iter++) { // loops until q([p],[p]')<=eps_q
		p = pi;
		//--------------------------------------------------------------------------
#if DIRECT == 1
		bool conv = kolevDirect(M1, b1, pi, w, t1, t3, t2, x0);
#else
		bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter); // outer p-solution
#endif
		//--------------------------------------------------------------------------
		// NOTE: we have to take innl and innu since the inner solution can be empty
		kolevInner(x0, t1, t3, t2, innl, innu); // inner solution to get bounds on the hull
		//--------------------------------------------------------------------------
		//xdm = w[m].inf(); xgm = w[m].sup(); // xdm corrsp. to xdi0, xgm corrsp. xgi0
		maxr = maxrad(pi);
		if (maxr < EPSP) {
			std::cout << "SUCCESS!!! in " << iter << " iterations" << std::endl;

			//% computing the final enclosure
#if DIRECT == 1
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
#else
			bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter); // outer p-solution
#endif
			//% computing the points at which the extrema are attained
			ivector pmid(K);
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk] << " ";
				pmid[kk] = pi[kk].mid();
			}
			std::cout << std::endl;
			if (minmax == 1) {
				std::cout << "max: " << w << std::endl;
			}
			else {
				std::cout << "min: " << w << std::endl;
			}
			//% computing the final solution
			eIntPLP(minmax, c, x0, t1, t2, t3, mw);
			std::cout << "eI = " << mw << std::endl;
			/*if (mws >= eI.inf() && mws <= eI.sup()) {
			std::cout << "RIGHT!!!" << std::endl;
			}
			else {
			std::cout << "WRONG!!!: error = " << ISLIB_ABS(mws - eI.mid()) << std::endl;
			}*/
			time = difftime(clock(), start) / CLOCKS_PER_SEC;
			return;
		}
		//%
		//% Contracting the interval using simple constraint propagation method
		contractorPLP(red, minmax, M1, b1, c, pi, innl, innu, t1, t3, t2, x0);
		//%
		if (red < EPSQ) {
			//%
			//% No improvement: only crude enclosure found
			//% HERE we can try to check the monotonicity 
			std::cout << "CRUDE!!! in " << iter << " iterations" << std::endl;
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].inf() << " ";
			}
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].sup() << " ";
			}
			std::cout << std::endl;
			/*bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			eIntPLP(minmax, c, x0, t1, t2, t3, mw);*/

			//std::cout << "pi" << std::endl << pi << std::endl;

			//% trying monotonicity
			MonotPLP(minmax, c, pi, M1, b1, mw);
			//std::cout << "CRUDE mmw(monot)" << mw << std::endl;
			return;
		}
	}
	std::cout << "FAILURE: maximum iterations exceeded!!!" << std::endl;
	return;
}

//%//////////////////////////
//% old version of contractor
//%--------------------------
void
contractorPLP2(real& red, int minmax, const dvector& c, ivector& p1, const dvector& innl, const dvector& innu, const dvector& t1, const dvector& t3,
const dmatrix& t2, const dvector& x0)
// Kolev's version slightly modified
{
	int i, n = t1.size(), K = p1.size(); // we optimize the m-th entry of the solution
	real d, t2sum = 0.0;
	dvector at2m(K), t2m(K);
	ivector p(K);

	omatrix oA(n, n);
	ovector ob(n);
	ivector q(1); q[0] = 0.0;

	p.fill_in(interval(-1.0, 1.0));

	for (i = 0; i < K; i++) {
		t2m[i] = 0.0;
		for (int m = 0; m < n; m++) {
			t2m[i] += c[m] * t2(m, i); // Lj
		}
		t2m[i] = t2m[i];
		at2m[i] = ISLIB_ABS(t2m[i]); // |Lj|
	}
	real t3m = 0.0;
	for (int m = 0; m < n; m++){
		t3m = t3m + ISLIB_ABS(c[m]) * t3[m]; // s
	}

#if (EVOLPLP==2)
	real sevo = 0.0, x0c = 0.0;
	for (int m = 0; m < n; m++) {
		x0c = x0c + c[m] * (x0[m] + t1[m]);
		//evolution(minmax, m, p1, q, oA, ob, inf); // trzeba sprawdzic, czy w ewolucujnym jest wszystko dobrze zdefiniowane
		//DE_PF(minmax, i, p1, q, oA, ob, inf);
		//sevo = sevo + minmax * c[m] * inf;
	}
	d = sevo - x0c; //(x0[m] + t1[m]);// tu trzeba odejmowac somu
#else
	if (minmax == -1) {
		d = t3m;
		for (int m = 0; m < n; m++) {
			d = d - at2m[m];
		}
	}
	else {
		d = -t3m;
		for (int m = 0; m < n; m++) {
			d = d + at2m[m];
		}
	}
#endif
	int kmax = 0;
	for (i = 1; i < K; i++) if (at2m[i] > at2m[kmax]) kmax = i;

	interval sum_t2p = 0.0;
	for (i = 0; i < K; i++) {
		if (i != kmax) {
			sum_t2p = sum_t2p + t2m[i] * p[i]; // c*L*p, t2 corresponds to L in the paper
		}
	}
	//-------------------------------------------------------------------------------

	interval bb = d - sum_t2p - interval(-t3m, t3m); // u Koleva po prostu [b]
	interval qq = bb / t2m[kmax];
	//-------------------------------------------------------------------------------

	interval pk = interval(ISLIB_MAX(qq.inf(), p[kmax].inf()), ISLIB_MIN(qq.sup(), p[kmax].sup()));
	real bkd = bb.inf();
	real bkg = bb.sup();
	real pkd = pk.inf();
	real pkg = pk.sup();
	real rpk = pk.rad(); // rk

	if ((pkd>-1.0) && (pkg < 1.0)) {
		// pKI strictly included in p1
		red = 1 - rpk;
		p1[kmax] = interval(p1[kmax].inf() + p1[kmax].rad()*(pkd + 1.0), p1[kmax].inf() + p1[kmax].rad()*(pkg + 1.0));;
		return;
	}
	else if ((pkd > 1.0) || (pkg < -1.0)) {
		// no intersection
		std::cout << "No intersection" << std::endl;
		red = -200;
		return;
	}
	else {
		real rpmax = p1[kmax].rad(); // r_k
		real rs = rpmax * rpk;
		real pc;
		red = 1.0 - rpk;
		if (pkd > -1.0) { // left reduction
			pc = p1[kmax].sup() - rs;
		}
		else if (pkg < 1.0) { // right reduction
			pc = p1[kmax].inf() + rs;
		}
		else {
			// no reduction
			red = 0;
			return;
		}
		p1[kmax] = interval(pc - rs, pc + rs);
	}
}
