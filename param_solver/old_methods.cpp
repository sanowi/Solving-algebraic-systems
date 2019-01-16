/**--------------------------------------------------*
 *  Different versions of Skalna's Direct Method
 *---------------------------------------------------*/
#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../interval/interval_base.h"
#include "../affine/aafr.h"
#include "../utils/randnum.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/direct_solvers.h"
#include "../truss/TrussStructure.h"
#include "old_methods.h"

static const real PGSIeps = 1.0e-8;
static const real eps = 1.0e-9;
static const real mi = spn(); //% smallest positive number

#define MAXITERMPO	100
#define MAXITER		200
#define EPSX		1.0e-11 //% accuracy

int
//%----------------------------------------------------------------------------
//% Checks if a revised affine matrix A is an H-matrix. If the method returns  
//% 'true', this means that A is an H-matrix, otherwise it is not decided.
//% Def: A is an H-matrix if there exists a real vector u > 0, such that 
//% <A>u > 0, where <A> - Ostrovsky comparison matrix the matrix A in argument 
//% is usually the Ostrovsky matrix
//%----------------------------------------------------------------------------
HMatrix(const dmatrix& A) //% revised affine matrix
{
	if (!A.is_square()) return NonSquareMatrix;

	int		n = A.num_rows();
	dvector	e(n), x(n);

	e.fill_in(1);
	if (gw_gausse(A, e, x)) { //% rozwiazanie ukladu <A>x = e
		for (int i = 0; i < n; i++) {
			if (!(x[i] > 0.0)) {
				return NotDecided;
			}
		}
		return Hmatrix; //% is H-matrix
	}
	return SingularMatrix; //% not H-matrix (singularity)
} //% >>> Hmatrix <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% with pre-conditioning using the midpoint inverse
//% The iterative scheme of the method is as follows:
//% x=D^{-1}(b-(L+U)x)
//% It origines from the following decompsition of the system matrix:
//% Ax=b <=> (L+D+U)x=b <=> Dx=b-(L+U)x <=> x=D^{-1}(b-(L+U)x)
//% http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
//%---------------------------------------------------------------------------- 
ParamGaussSeidelNV(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector w1(n), c(n);
	dvector xinf(n), xsup(n);
	ivector xp(n);
	aafrmatrix C(n, n);
	dmatrix R(n, n);

	//% initial enclosure of the pre-conditioned system
	err = ParametricDirectMethod(A, b, xp, R); //% could be changed to Rump's method
	if (err != Success) return err;
	//% transform initial interval enclosure into revised affine form
	//% with midpoint of interval as central value and radius of interval
	//% as radius of accumulative error, no noise variable is introduced
	w.fill_in(0.0);
	for (int i = 0; i < n; i++) { //% O(n)
		w[i] = w[i] + xp[i];
	}
	//% left pre-conditioning
	//xp.fill_in(0.0);
	C = R * A; //% O(m*n^3*K) (m-cost of multiplication of revised affine forms)
	c = R * b; //% O(n^2*K)
			   //%-------------------------------------------------------------
			   //% Main iteration loop
			   //%-------------------------------------------------------------
	niter = 0;
	do {
		niter++;
		if (niter > MAXITER) break;
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
	} while (dist1(w, w1) > EPSX);
	return Success;
} //% >>> ParamGaussSeidel <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Gauss-Seidel iteration <<<
//% with right pre-conditioning with midpoint inverse
//%----------------------------------------------------------------------------
ParamGaussSeidelIINV(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	aafrvector& w, //% p-solution (revised affine vector)
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector w1(n);
	dvector xinf(n), xsup(n);
	ivector xp(n);
	aafrmatrix B(n, n);
	dmatrix	R(n, n), Am(n, n);

	err = ParametricDirectMethodPII(A, b, xp, R); //% initial solution
	if (err != Success) return err;
	//% transform initial interval enclosure into revised affine form
	//% with midpoint of interval as central value and radius of interval
	//% as radius of accumulative error, no noise variable is introduced
	w.fill_in(0.0);
	for (int i = 0; i < n; i++) { // O(n)
		w[i] = w[i] + xp[i];
	}
	B = A * R; //% right pre-conditioning
	niter = 0;
	//%-------------------------------------------------------------
	//% Main iteration loop
	//%-------------------------------------------------------------
	do {
		niter++;
		if (niter > MAXITER) { return NotHMatrix; }
		w1 = w;
		for (int i = 0; i < n; i++) { // all loop O(n^2*K*log(K)), O(n^2*m), m-cost of multiplication
			AAFR s = 0.0;
			for (int j = 0; j < i; j++) { //% lower triangular
				if (B(i, j).reduce() != 0.0 && w[j].reduce() != 0.0) {
					s = s + B(i, j) * w[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			for (int j = i + 1; j < n; j++) { //% upper triangular
				if (B(i, j).reduce() != 0.0 && w1[j].reduce() != 0.0) {
					s = s + B(i, j) * w1[j]; //% O(K*log(K)), O(m), m-cost of multiplication
				}
			}
			w[i] = (b[i] - s) / B(i, i); //% O(K*log(K))
		}
	} while (dist1(w, w1) > EPSX);
	w = R * w; //% final solution
	return true;
} //% >>> ParamGaussSeidelII <<<

bool
//%----------------------------------------------------------------------------
//% Skalna's Method for solving PILS
//% with one vector of interval parameters
//% - p - vector of left hand interval parameters
//% - q - vector of right hand interval parameters
//% - oA - matrix of left hand dependencies (A = oA * p)
//% - ob - vector of right hand dependencies (b = ob * p)
//%----------------------------------------------------------------------------
ISM(const ivector& p, 
	const omatrix& oA, 
	const ovector& ob, 
	ivector& w)
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
	dvector		x(n), bm(n), cmz(n), magz(n), pm(np);
	dmatrix		R(n, n), Cm(n, n);
	ivector		z(n);
	imatrix		C(n, n);

	w = ivector(n);
	mid(p, pm);
	bm = VectorFromDep(pm, ob);
	R = MatrixFromDep(pm, oA);
	if (inv(R)) {
		x = R * bm;
		c_matrix(R, p, oA, C);
		mig(C, Cm);
		if (!HMatrix(Cm)) { //% H-matrix property verification
			return false; //% not an H-matrix
		}
		z_vector(x, R, p, oA, ob, z); //% z = z(p) = R(b(p) - A(p)x0)
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0); //interval(-cmz[i], cmz[i]);
			}
			return true;
		}
	}
	std::cout << "Singular matrix!!!" << std::endl;
	return false; // singular matrix
}

bool
//%----------------------------------------------------------------------------
//% Skalna's Method for solving PILS
//% with two vectors of interval parameters
//% - p - vector of left hand interval parameters
//% - q - vector of right hand interval parameters
//% - oA - matrix of left hand dependencies (A = oA * p)
//% - ob - vector of right hand dependencies (b = ob * p)
//%----------------------------------------------------------------------------
ISM(const ivector& p,
	const ivector& q,
	const omatrix& oA,
	const ovector& ob,
	ivector& w)
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq), pm(np);
	dmatrix		R(n, n), Cm(n, n);
	ivector		z(n);
	imatrix		C(n, n);

	//%https://msdn.microsoft.com/en-us/library/vstudio/dd492418%28v=vs.100%29.aspx
	w = ivector(n);
	mid(q, qm);
	mid(p, pm);
	
	bm = VectorFromDep(qm, ob);
	R = MatrixFromDep(pm, oA);
	if (inv(R)) {
		x = R * bm;
		c_matrix(R, p, oA, C);
		mig(C, Cm);
		if (!HMatrix(Cm)) { //% H-matrix property verification
			return false; //% not an H-matrix
		}
		z_vector(x, R, p, q, oA, ob, z); //% z = z(p, q) = R(b(q) - A(p)x0)
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0);
			}
			return true;
		}
	}
	return false; //% singular matrix
}

bool
ISM(const ivector& p, const ivector& q, const dmatrix& R, const imatrix& C, const omatrix& oA,
	const ovector& ob, ivector& w)
	//%----------------------------------------------------------
	//% ISM version with two sets of interval vectors
	//% the difference is that R=(mid(A))^{-1} and C = <RA>
	//% are passed as arguments
	//%-----------------------------------------------------------
{
	int			i, n = oA.num_rows(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq);
	dmatrix		Cm(n, n);
	ivector		z(n);

	w = ivector(n);
	mid(q, qm);
	bm = VectorFromDep(qm, ob);
	x = R * bm;
	mig(C, Cm);
	if (!HMatrix(Cm)) { //% H-matrix property verification
		return false; //% not an H-matrix
	}
	z_vector(x, R, p, q, oA, ob, z);
	mag(z, magz);
	if (gw_gausse(Cm, magz, cmz)) {
		for (i = 0; i < n; i++) {
			w[i] = x[i] + cmz[i] * interval(-1.0, 1.0); //interval(-cmz[i], cmz[i]); //;
		}
		return true;
	}
	return false; //% singular matrix
}

bool
ISMC(const ivector& p, const ivector& q, const dmatrix& R, const imatrix& C, const omatrix& oA,
	const ovector& ob, ivector& w)
	//%--------------------------------------------------------
	//% ISMF version with C = <RA> and R are being passed 
	//% as an argument and with computing b = ob * q
	//%--------------------------------------------------------
{
	if (!oA.is_square()) return false;
	int			i, n = oA.num_rows(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq);
	dmatrix		Cm(n, n);
	ivector		b(n), z(n);

	w = ivector(n);
	mid(q, qm);
	bm = VectorFromDep(qm, ob); // tu zmiana 01.03.2009 - z pv na qv
	x = R * bm;
	mig(C, Cm);
	if (!HMatrix(Cm)) { // H-matrix property verification
		return false; // not an H-matrix
	}
	z_vector(x, R, p, q, oA, ob, z);
	mag(z, magz);
	if (gw_gausse(Cm, magz, cmz)) {
		for (i = 0; i < n; i++) {
			w[i] = x[i] + cmz[i] * interval(-1, 1);
		}
		return true;
	}
	return false; //% singular matrix
}

int
//%----------------------------------------------------------------------------
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The method is based on Krawczyk iteration and computes verified
//% outer bounds for the united parametric solution set.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	ivector& w, //% interval solution vector
	ivector& ww, //% inner estimation of hull solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), x0(n);
	dmatrix	Ac(n, n), Cm(n, n), R(n, n), I(n, n);
	imatrix	B(n, n), C(n, n);
	ivector	z(n), y0(n), y1(n), yk(n), yy(n), imi(n), d(n);
	interval ykk;

	mid(A, Ac); //% midpoint matrix
	R = Ac;
	if (inv(R)) { //% approximate midpoint inverse
		mid(b, bc); //% midpoint vector
		x0 = R * bc; //% midpoint solution
		//ResidIter(R, Ac, bc, x0);
		Unit(I);
		B = I - reduce(R*A); //% [B] = I - [R*A(p)]
		mig(B, Cm);
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		z = reduce(R*(b - A * x0));
		y0 = z; //% starting point for the iteration
		niter = 0;
		//% Main iteration loop (verification step)
		while (true) {
			niter++;
			if (niter > MAXITERMPO) return Failed;
			blow(y0, eps, mi);
			y1 = y0;
			for (int i = 0; i < n; i++) {
				ykk = 0.0;
				for (int j = 0; j < n; j++) {
					ykk = ykk + B(i, j) * y1[j];
				}
				y1[i] = z[i] + ykk;
			}
			if (y0.set_strictly_contains(y1)) { //% solution found
				d = B * y1;
				w = x0 + y1;
				//% inner estimation of hull
				for (int i = 0; i < n; i++) {
					ww[i] = x0[i] + interval(z[i].inf() + d[i].sup(), z[i].sup() + d[i].inf());
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
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The method uses specific representation of parametric systems,
//% i.e., representation A(p)=\sum_{k=1}^KA^k\cdot p_k
//%----------------------------------------------------------------------------
Rump(const ivector& p, //% vector of interval parameters
	const omatrix& oA, //% matrix of left hand dependencies
	const ovector& ob, //% vector of right hand dependencies
	ivector& w, //% interval solution vector
	int& stepCount) //% number of iterations
{
	if (!oA.is_square()) return NonSquareMatrix;

	if (p.size() != oA(0, 0).size() || oA.num_rows() != ob.size()) return false;
	int	n = oA.num_rows(), np = p.size();
	dvector	bm(n), x0(n), yylb(n), yyub(n), pm(np), y0r(n);
	dmatrix	Cm(n, n), R(n, n), I(n, n);
	imatrix	A(n, n), B(n, n), C(n, n);
	ivector	b(n), z(n), y0(n), y1(n), yk(n), yy(n), imi(n);

	mid(p, pm);
	R = MatrixFromDep(pm, oA); //% midpoint matrix
	bm = VectorFromDep(pm, ob); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		x0 = R * bm; //% midpoint solution
		c_matrix(R, p, oA, C); //% C = R*A(p)
		Unit(I);
		B = I - C;
		mig(B, Cm); // Cm = <C>
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		z_vector(x0, R, p, oA, ob, z); //% z = R*(b(p)-A(p)*x0)
		y0 = z;
		stepCount = 0;
		//% Main iteration loop
		while (true) {
			stepCount++;
			if (stepCount > MAXITERMPO) return Failed;
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
//% Self-verified fixed point iteration for solving parametric interval 
//% systems. The computation is performed by using bpth intervals and
//% revised affine forms.
//%----------------------------------------------------------------------------
RumpAAFR(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised right-hand side vector
	ivector& w, //% interval outer solution
	ivector& ww, //% interval inner solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	aafrvector z(n), y1(n), yk(n), imi(n), d(n);
	aafrmatrix B(n, n), C(n, n); //% this is the main difference
	cvector xpidx(1);
	dvector	bm(n), x0(n), xpcfs(1);
	ivector y0(n), y1i(n), di(n), zi(n);
	dmatrix	Cm(n, n), R(n, n), I(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		mid(b, bm); //% midpoint vector
		x0 = R * bm; //% midpoint solution
		C = R * A;
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = I(i, j) - C(i, j);
			}
		}
		mig(reduce(B), Cm);
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		z = R * (b - A * x0);
		y0 = reduce(z); //% initial solution
		niter = 0;
		//% Main iteration loop
		while (true) {
			niter++;
			blow(y0, eps, mi);
			xpidx[0] = 0;
			for (int i = 0; i < n; i++) { //% transform: interval => aafr
				xpcfs[0] = y0[i].mid();
				y1[i] = AAFR(xpidx, xpcfs, y0[i].rad());
			}
			for (int i = 0; i < n; i++) {
				AAFR ykk = 0.0;
				for (int j = 0; j < n; j++) {
					ykk = ykk + B(i, j) * y1[j];
				}
				y1[i] = z[i] + ykk;
			}
			y1i = reduce(y1);
			if (y0.set_strictly_contains(y1i)) {
				d = B * y1;
				di = reduce(d);
				w = x0 + reduce(y1);
				zi = reduce(z);
				//% inner estimation of hull solution
				for (int i = 0; i < n; i++) {
					ww[i] = x0[i] + interval(zi[i].inf() + di[i].sup(), zi[i].sup() + di[i].inf());
				}
				return Success;
			}
			y0 = reduce(y1);
		}
	}
	return SingularMatrix;
} //% >>> RumpAAFR <<<

int
//%----------------------------------------------------------------------------
//% >>> Rump for systems with Multiple Right-hand Side (MHRS) <<<
//% Self-verified fixed point iteration for solving parametric interval 
//% systems with multiple righ-hand side. The method is based on Krawczyk 
//% iteration.
//%----------------------------------------------------------------------------
Rump(const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& Bb, //% revised affine right-hand side matrix
	imatrix& w, //% interval solution vector
	imatrix& ww, //% inner estimation of hull solution
	int& niter) //% number of iterations
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), m = Bb.num_cols();
	dmatrix	Bc(n, m), X0(n, m);
	dmatrix	Cm(n, n), R(n, n), I(n, n);
	imatrix	B(n, n), C(n, n), imi(n, m);
	imatrix	Z(n, m), Y0(n, m), Y1(n, m), yk(n, m), yy(n, m), D(n, m);
	aafrmatrix Ax(n, m);
	interval ykk;

	mid(A, R);
	if (inv(R)) {
		mid(Bb, Bc);
		X0 = R * Bc;
		C = reduce(R*A);
		Unit(I);
		B = I - C;
		mig(B, Cm);
		if (rhoSpectral(Cm) >= 1.0) return NotStronglyRegular;
		Z = reduce(R*(Bb - A * X0));
		Y0 = Z;
		niter = 0;
		while (true) {
			niter++;
			blow(Y0, eps, mi);
			Y1 = Y0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					ykk = 0.0;
					for (int k = 0; k < n; k++) {
						ykk = ykk + B(i, k) * Y1(k, j);
					}
					Y1(i, j) = Z(i, j) + ykk;
				}
			}
			if (StrictlyContains(Y1, Y0)) {
				D = B * Y1;
				w = X0 + Y1;
				//% inner estimation of hulll
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						ww(i, j) = X0(i, j) + interval(Z(i, j).inf() +
							D(i, j).sup(), Z(i, j).sup() + D(i, j).inf());
					}
				}
				return Success;
			}
			Y0 = Y1;
		}
	}
	return SingularMatrix;
} //% >>> Rump for MHRS <<<

static dvector
FromDepend(const dvector& p, const ovector& ob)
{
	int n = ob.size(), K = p.size();
	dvector b(n);

	b = ob[0];
	for (int k = 1; k < K + 1; k++) {
		b = b + ob[k] * p[k - 1];
	}
	return b;
}

//%==============================================================
//% procedures from "zvector" file
//%==============================================================

bool
c_matrix(const dmatrix& R, const ivector& p, const omatrix& oA, imatrix& C)
//% Calculates [C] = R * A([p])
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
	dvector		s(np);
	//-----------------------------
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			s.fill_in(0);
			for (int k = 0; k < n; k++) {
				s += oA(k, j) * R(i, k);
			}
			C(i, j) = p * s;
		}
	}
	return true;
}

bool
c_matrix(const imatrix& R, const ivector& p, const omatrix& oA, imatrix& C)
//% Calculates [C] = [R] * A(p)
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
	ivector		s(np);
	//-----------------------------
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			s.fill_in(0);
			for (int k = 0; k < n; k++) {
				s += oA(k, j) * R(i, k);
			}
			C(i, j) = p * s;
		}
	}
	return true;
}

bool
c_matrix(const omatrix& B, const ivector& p, imatrix& C)
//% Calculates [C] based on [B] = R * oA
{
	if (!B.is_square()) return false;
	int n = B.num_rows();
	//-----------------------------
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			C(i, j) = p * B(i, j);
		}
	}
	return true;
}

bool
b_matrix(const dmatrix& R, const ivector& p, const omatrix& oA, omatrix& B)
//% Calculates [B]  = R * oA 
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
	dvector		s(np);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			s.fill_in(0);
			for (int k = 0; k < n; k++) {
				s += oA(k, j) * R(i, k);
			}
			B(i, j) = s;
		}
	}
	return true;
}

bool
b_vector(const dmatrix& R, const ivector& q, const ovector& ob, ivector& b)
//% Calculates [C] = R * A([p])
{
	if (R.num_cols() != ob.size()) return false;
	int			n = ob.size();
	dvector		s(q.size());

	for (int i = 0; i < n; i++) {
		s.fill_in(0);
		for (int k = 0; k < n; k++) {
			s += ob[k] * R(i, k);
		}
		if (s != 0.0) {
			b[i] = q * s;
		}
		else {
			b[i] = 0;
		}
	}
	return true;
}

bool
z_vector(const dvector& x, const dmatrix& R, const ivector& p, const ivector& q, const omatrix& oA,
	const ovector& ob, ivector& z)
	//%--------------------------------------------
	//% Calculates [z] = R*([b(q)] - [A(p)]*x0)
	//% R - real matrix, inverse of midpoint
	//% p - left hand Vector of parameters
	//% q - right hand Vector of parameters
	//%--------------------------------------------
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size(), nq = q.size();
	dvector		s1(nq), s2(np);

	for (int i = 0; i < n; i++) {
		s1.fill_in(0);
		s2.fill_in(0);
		for (int j = 0; j < n; j++) {
			s1 += ob[j] * R(i, j);
			for (int k = 0; k < n; k++) {
				s2 += oA(j, k) * (R(i, j) * x[k]);
			}
		}
		z[i] = (q * s1) - (p * s2); // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
	}
	return true;
}

bool
z_vector(const ivector& x, const imatrix& R, const ivector& p, const omatrix& oA, const ovector& ob, ivector& z)
//%---------------------------------------------
//% Calculates [z] = R*([b(q)] - [A(p)]*x0)
//% R - interval matrix, inverse of midpoint
//% p - left hand Vector of parameters
//% q - right hand Vector of parameters
//%---------------------------------------------
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size();
	ivector		s1(np), s2(np);

	for (int i = 0; i < n; i++) {
		s1.fill_in(0.0);
		s2.fill_in(0.0);
		for (int j = 0; j < n; j++) {
			s1 = s1 + ob[j] * R(i, j);
			for (int k = 0; k < n; k++) {
				s2 = s2 + oA(j, k) * (R(i, j) * x[k]);
			}
		}
		z[i] = p * (s1 - s2);
	}
	return true;
}

bool
z_vector(const ivector& x, const imatrix& R, const ivector& p, const ivector& q, const omatrix& oA,
	const ovector& ob, ivector& z)
	//%---------------------------------------------
	//% Calculates [z] = R*([b(q)] - [A(p)]*x0)
	//% R - interval matrix, inverse of midpoint
	//% p - left hand Vector of parameters
	//% q - right hand Vector of parameters
	//%---------------------------------------------
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size(), nq = q.size();
	ivector		s1(nq), s2(np);
	//-----------------------------
	for (int i = 0; i < n; i++) {
		s1.fill_in(0.0);
		s2.fill_in(0.0);
		for (int j = 0; j < n; j++) {
			s1 = s1 + ob[j] * R(i, j);
			for (int k = 0; k < n; k++) {
				s2 = s2 + oA(j, k) * (R(i, j) * x[k]);
			}
		}
		z[i] = q * s1 - p * s2; // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
	}
	return true;
}

bool
z_vector(const dvector& b, const dvector& x, const dmatrix& R, const ivector& p, const omatrix& oA,
	const ovector& ob, ivector& z)
	//%---------------------------------------------
	//% Calculates [z] = R*([b(p)] - [A(p)]*x0)
	//% p - left hand Vector of parameters
	//%---------------------------------------------
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
	dvector		s1(np), s2(np);

	for (int i = 0; i < n; i++) {
		s1.fill_in(0);
		s2.fill_in(0);
		for (int j = 0; j < n; j++) {
			s1 += ob[j] * R(i, j);
			for (int k = 0; k < n; k++) {
				s2 += oA(j, k) * x[k] * R(i, j);
			}
		}
		z[i] = p * (s1 - s2);
	}
	return true;
}

bool
z_vector(const dvector& x, const dmatrix& R, const ivector& p, const omatrix& oA, const ovector& ob, ivector& w)
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size();
	ivector		s1(np), s2(np);

	for (int i = 0; i < n; i++) {
		s1.fill_in(0);
		s2.fill_in(0);
		for (int j = 0; j < n; j++) {
			s1 = s1 + ob[j] * R(i, j);
			for (int k = 0; k < n; k++) {
				s2 = s2 + oA(j, k) * (R(i, j) * x[k]);
			}
		}
		w[i] = p * (s1 - s2);
	}
	return true;
}

//===============================================================================================
// RECTANGULAR SYSTEMS
//===============================================================================================

bool ISMO(const ivector& p, const omatrix& oA, const ovector& ob, ivector& w)
//%-----------------------------------------------------------------
//% Solving OVER- and UNDERdetermined systems using Rump's method
//% Skalna's Method for solving PILS
//% with one set of interval parameters
//% - p - vector of left hand interval parameters
//% - oA - matrix of left hand dependencies (A = oA * p)
//% - ob - vector of right hand dependencies (b = ob * p)
//%-----------------------------------------------------------------
{
	int	n = oA.num_rows(), m = oA.num_cols(), np = p.size(), nm;
	if (n == m) {
		nm = n;
	}
	else {
		nm = n + m;
	}
	dvector		x(nm), bm(nm), cmz(nm), magz(nm), pm(np), el(np);
	dmatrix		R(nm, nm), Cm(nm, nm);
	ivector		z(nm);
	imatrix		C(nm, nm);
	omatrix		oAA(nm, nm);
	ovector		obb(nm);
	w = ivector(m);
	if (n == m) {
		oAA = oA;
	}
	else {
		el.fill_in(0.0);
		oAA.fill_in(el);
		obb.fill_in(el);
		if (n > m) { // overdetermined
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					oAA(i, j) = oA(i, j);
				}
				el[0] = -1.0;
				oAA(i, i + m) = el;
				obb[i] = ob[i];
			}
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					oAA(i + n, j + m) = oA(j, i);
				}
			}
		}
		else { // underdetermined
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					oAA(i, j) = oA(j, i);
				}
				el[0] = -1.0;
				oAA(i, i + n) = el;
				obb[i + m] = ob[i];
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					oAA(i + m, j + n) = oA(i, j);
				}
			}
		}
	}

	mid(p, pm);
	bm = VectorFromDep(pm, obb);
	R = MatrixFromDep(pm, oAA);
	if (inv(R)) {
		x = R * bm;
		c_matrix(R, p, oAA, C);
		mig(C, Cm);
		if (!HMatrix(Cm)) { // H-matrix property verification
			std::cout << "Probably not an H-matrix" << std::endl;
			return false; // not an H-matrix
		}
		z_vector(x, R, p, oAA, obb, z); /* z = z(p) = R(b(p) - A(p)x0) */
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < m; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0); //interval(-cmz[i], cmz[i]);
			}
			return true;
		}
	}
	std::cout << "Singular matrix!!!" << std::endl;
	return false; // singular matrix
}