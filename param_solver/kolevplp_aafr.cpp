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
#include "../examples/pexamples.h"
#include "kolev.h"
#include "kolevplp_aafr.h"

#define EPSP 1.0e-3
#define EPSQ 1.0e-3
#define MAXITER 1000

#define EVOLPLP 0
#define tau1 0.0 //right // for rho=0.1,...,0.263 tau1=0.18, tau2=0.02
#define tau2 0.0 //left
#define tau11 0.5 //4 //0.9
#define tau22 0.5 //02 // 0.02

void 
kolevPLPOuter(const dvector& c, const aafrmatrix& M1, const aafrvector& b1, interval& mw, real& time)
// calculates m-th entry of the hull solution
{
	int n = b1.size(), niter, precond = 0;
	interval eI;
	aafrvector w(n);
	real x0 = 0.0, s = 0.0, xr = 0.0, xi = 0.0;

	//--------------------------------------------------------------------------
	bool isok = ParamGaussSeidel(M1, b1, w, niter); // outer p-solution
	if(isok) std::cout << "OK" << std::endl;
	/*for(int i = 0; i < n; i++) {
		std::cout << w[i].reduce() << std::endl;
	}*/
	int np = w[0].cfs().size();
	for(int j = 1; j < np; j++) {
		s = 0.0;
		for(int i = 0; i < n; i++) {
			s += w[i].cfs()[j];
		}
		xi += ISLIB_ABS(s);
	}
	for(int i = 0; i < n; i++) {
		x0 += w[i].cfs()[0];
		xr += w[i].r();
	}
	mw = x0 + interval(-xi-xr, xi+xr);
	return;
}

void
contractorPLP(real& red, int minmax, aafrmatrix& M, aafrvector& b, const dvector& c, ivector& p1, 
const dvector& innl, const dvector& innu, const dvector& t1, const dvector& t3, const dmatrix& t2, const dvector& x0)
//%
//% Kolev's version slightly modified
//%
{
	int i, n = t1.size(), K = p1.size(), niter, precond = 0; // we optimize the m-th entry of the solution
	real t2sum = 0.0;
	interval d; // d is like l^u in the paper
	dvector at2m(K), t2m(K);
	ivector p(K);

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

	//--------------------
	// Example 3
	//real tauL = 0.74;
	//real tauU = 0.71; // 0.71 is good for 0.3, but not good for 0.4

	real tauL = 0.8;
	real tauU = 0.94;
	ivector pt(K), ww(n);
	aafrvector wa(n);
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
	aafrvector ptaafr(pt.size());
	for(int i = 0; i < pt.size(); i++) {
		ptaafr[i] = AAFR(pt[i]);
	}
	//%-----------------
	//% trzeba stworzyc nowa macierz M i nowy wektor b tak jak w przykladzie odpowiednim
	//%-----------------
	Ex3(M, b, ptaafr, K, 0.1);
	//%-----------------
	ParamGaussSeidel(M, b, wa, niter);
	ww = reduce(wa);
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

	//--------------------

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
		p1[kmax] = interval(p1[kmax].inf() + p1[kmax].rad()*(pkd + 1.0), p1[kmax].inf() + p1[kmax].rad()*(pkg + 1.0));
		std::cout << "reduction1: " << p1[kmax] << std::endl;
		char c = _getch();
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
		std::cout << "reduction2: " << p1[kmax] << std::endl;
		char c = _getch();
	}
}

#define DIRECT 0

void 
kolevPLP(int minmax, const dvector& c, const aafrmatrix& M1, const aafrvector& b1, interval& mw, int& iter, real& time)
//%
//% calculates m-th entry of the hull solution
//% M1 and b1 are generated in testMKM from the selected example
{
	int n = b1.size(), K = b1[0].idx().size(), niter, precond = 0;
	interval eI;
	real maxr, red;

	dvector t1(n), t3(n), pp(K), bi(n), xa(n), innl(n), innu(n), b0(n), x0(n);
	dmatrix t2(n, K), A0(n, n), A(n, n);

	aafrvector w(n);
	ivector p(K), pi(K), hx(n), iw(n), ww(n);
	aafrmatrix M(n,n);
	aafrvector b(n);

	time_t start = clock();

	M = M1;
	b = b1;

	for (iter = 0; iter < MAXITER; iter++) { // loops until q([p],[p]')<=eps_q
		pi.fill_in(interval(-1.0, 1.0));
		//--------------------------------------------------------------------------
		bool conv = ParamGaussSeidel(M, b, w, niter); // outer p-solution
		//--------------------------------------------------------------------------
		// NOTE: we have to take innl and innu since the inner solution can be empty
		IterInner(w, innl, innu); // inner solution to get bounds on the hull
		//--------------------------------------------------------------------------
		//xdm = w[m].inf(); xgm = w[m].sup(); // xdm corrsp. to xdi0, xgm corrsp. xgi0
		iw = reduce(w);
		std::cout << "wouter = " << std::endl << ww << std::endl;
		maxr = maxrad(iw);
		if (maxr < EPSP) {
			std::cout << "SUCCESS!!! in " << iter << " iterations" << std::endl;

			//% computing the final enclosure

			bool conv = ParamGaussSeidel(M, b, w, niter); // outer p-solution

			//% computing the points at which the extrema are attained
			/*for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].inf() << " ";
			}
			for (int kk = 0; kk < K; kk++) {
				std::cout << pi[kk].sup() << " ";
			}
			std::cout << std::endl;*/
			
			//% computing the final solution
			real x0 = 0.0, s = 0.0, xr = 0.0, xi = 0.0;
			int np = w[0].cfs().size();
			for(int j = 1; j < np; j++) {
				s = 0.0;
				for(int i = 0; i < n; i++) {
					s += w[i].cfs()[j];
				}
				xi += ISLIB_ABS(s);
			}
			for(int i = 0; i < n; i++) {
				x0 += w[i].cfs()[0];
				xr += w[i].r();
			}
			mw = x0 + interval(-xi-xr, xi+xr);
			
			std::cout << "eI = " << mw << std::endl;
			time = difftime(clock(), start) / CLOCKS_PER_SEC;
			return;
		}
		//%
		//% Contracting the interval using simple constraint propagation method
		contractorPLP(red, minmax, M, b, c, pi, innl, innu, t1, t3, t2, x0);
		//%
		if (red < EPSQ) {
			//%
			//% No improvement: only crude enclosure found
			//% HERE we can try to check the monotonicity 
			std::cout << "CRUDE!!! in " << iter << " iterations" << std::endl;
			/*bool conv = kolevOuter(M1, b1, pi, w, t1, t3, t2, x0, niter);
			eIntPLP(minmax, c, x0, t1, t2, t3, mw);*/
			
			//std::cout << "pi" << std::endl << pi << std::endl;
			
			//% trying monotonicity
			//monotPLP(minmax, c, pi, M1, b1, mw);
			//std::cout << "CRUDE mmw(monot)" << mw << std::endl;
			return;
		}
	}
	std::cout << "FAILURE: maximum iterations exceeded!!!" << std::endl;
	return;
}
