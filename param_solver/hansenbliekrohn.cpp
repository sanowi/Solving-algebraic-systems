//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Hansen-Bliek-Rohn method
//%
//%##########################################################################
//%

#include "../utils/stdafx.h"
#include "../affine/aafr.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/direct_solvers.h"
#include "../param_solver/old_methods.h"
#include "../truss/TrussStructure.h"
#include "hansenbliekrohn.h"

//%============================================================================
//%
//% Methods for solving paramteric interval linear systems
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% The method is a straightforward extension of the HBR method for ILS
//%----------------------------------------------------------------------------
HansenBliekRohn1(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n);
	dmatrix	R(n, n), Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = R * bc; //% midpoint solution
		z = R * reduce(b); //% left pre-conditioning of interval matrix
		C = R * reduce(A); //% left pre-conditioning of interval vector
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(xc, xcm); //% abs of midpoint solution
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (xcm + zr); //% x^*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2); //% lower bound
				w[i] = interval(wdown[i], wup[i]); //% i-th element of the solution
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn1 <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% First, the system is preconditioned with the midpoint inverse.
//% Then, the formulas from Rohn's paper are applied.
//% The method uses precomputed midpoint inverse to decreas time complexity.  
//% The method is used mainly in HBR method with residual correction.
//%----------------------------------------------------------------------------
HansenBliekRohn1(const dmatrix& R, //% midpoint inverse
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), zr(n), wdown(n), wup(n), xcm(n), xc(n), xast(n);
	dmatrix	Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(b, bc); //% midpoint vector
	xc = R * bc; //% midpoint solution
	mag(xc, xcm); //% abs of midpoint solution
	C = R * reduce(A); //% preconditioning
	rad(C, Cr); //% radius of C
	if (rhoSpectral(Cr) >= 1.0) {
		return NotStronglyRegular;
	}
	z = R * reduce(b);
	rad(z, zr); //% radius of z
	Unit(I);
	M = I - Cr;
	if (inv(M)) { //% inverse of M
		xast = M * (xcm + zr); //% x*
		for (int i = 0; i < n; i++) {
			real ti = (1.0 / (2.0 * M(i, i) - 1.0));
			real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
			real v2 = ti * v1;
			wup[i] = ISLIB_MAX(v1, v2); //% upper bound
			v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
			v2 = ti * v1;
			wdown[i] = ISLIB_MIN(v1, v2); //% lower bound
			w[i] = interval(wdown[i], wup[i]); //% i-th element of the solution
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn1 <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<< 
//% First, residual correction is applied, then the system is preconditioned 
//% with the midpoint inverse. Next, the formulas from Rohn's paper are applied.
//%----------------------------------------------------------------------------
HansenBliekRohn1RC(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	aafrvector z(n);
	dvector	bc(n), xc(n);
	dmatrix	R(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = R * bc; //% midpoint solution
		z = b - A * xc; //% residual correction
		err = HansenBliekRohn1(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn1RC <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% First, the system is preconditioned with the midpoint inverse.
//% Then, the formulas from Rohn's paper are applied.
//%----------------------------------------------------------------------------
HansenBliekRohn(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), compare = 1;
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n), errl(n), erru(n);
	dmatrix	R(n, n), Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);
	
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = R * bc; //% midpoint solution
		z = reduce(R * b); //% left pre-conditioning
		C = reduce(R * A); //% left pre-conditioning
		rad(C, Cr); //% radius of C
		std::cout << rhoSpectral(Cr) << std::endl;
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(xc, xcm); //% abs of midpoint solution
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (xcm + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2); //% lower bound
				w[i] = interval(wdown[i], wup[i]); //5 i-th element of the solution
			}
			if (compare == 1) {
				BS_vs_HBR(M, xc, xast, errl, erru);
				for (int i = 0; i < n; i++) {
					std::cout << errl[i] << ",  " << erru[i] << std::endl;
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% First, the system is preconditioned with the midpoint inverse.
//% Then, the formulas from Rohn's paper are applied.
//% The method uses precomputed midpoint inverse to decreas time complexity.  
//% The method is used mainly in HBR method with residual correction.
//%----------------------------------------------------------------------------
HansenBliekRohn(const dmatrix& R, //% midpoint inverse
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), wdown(n), wup(n), xcm(n), xc(n), xast(n), zr(n);
	dmatrix	Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(b, bc); //% midpoint vector
	xc = R * bc; //% midpoint solution
	mag(xc, xcm); //% abs of midpoint solution
	C = reduce(R * A); //% radius of C
	rad(C, Cr); //% radius of C
	if (rhoSpectral(Cr) >= 1.0) {
		return NotStronglyRegular;
	}
	z = reduce(R * b);
	rad(z, zr); //% radius of z
	Unit(I);
	M = I - Cr;
	if (inv(M)) { //% inverse of M
		xast = M * (xcm + zr); //% x*
		for (int i = 0; i < n; i++) {
			real ti = (1.0 / (2.0 * M(i, i) - 1.0));
			real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
			real v2 = ti * v1;
			wup[i] = ISLIB_MAX(v1, v2); //% upper bound
			v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
			v2 = ti * v1;
			wdown[i] = ISLIB_MIN(v1, v2); //% lower bound
			w[i] = interval(wdown[i], wup[i]);
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<< 
//% First, residual correction is applied, then the system is preconditioned 
//% with the midpoint inverse. Next, the formulas from Rohn's paper are applied.
//%----------------------------------------------------------------------------
HansenBliekRohnRC(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	aafrvector z(n);
	dvector	bc(n), xc(n);
	dmatrix	R(n, n);

	mid(A, R); //% midpoint matrix
	mid(b, bc); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse
		xc = R * bc; //% midpoint solution
		z = b - A * xc; //% residual correction
		err = HansenBliekRohn(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnRC <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<< 
//% with preconditioning and optionally with or without residual correction
//%----------------------------------------------------------------------------
HansenBliekRohnP(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	int r, //% 0 - without rc, 1 - with rc
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	aafrvector z(n);
	aafrmatrix C(n, n);
	dvector	bc(n), xc(n);
	dmatrix	R(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse - preconditioner
		mid(b, bc); //% midpoint right-hand vector
		xc = R * bc; //% midpoint solution
		C = R * A; //% left pre-conditioning
		if (r == 0) {
			z = R * b; //% left pre-conditioning
			return HansenBliekRohn(C, z, w);
		}
		z = R * (b - A * xc);
		err = HansenBliekRohn(C, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix; //% singular matrix
} //% >>> HansenBliekRohnRC <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<< 
//% The method is used in Krawczyk iteration to obtain the initital solution.
//% The method returns the midpoint inverse to decrease time complexity
//% of Krawczyk iteration.
//%----------------------------------------------------------------------------
HansenBliekRohn(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w, //% interval vector solution
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), mxc(n), xast(n);
	dmatrix	Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	mid(b, bc); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse
		xc = R * bc; //% midpoint solution
		mag(xc, mxc); //% abs of midpoint solution
		C = reduce(R * A); //% left pre-conditioning
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		z = reduce(R * b);
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (mxc + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - mxc[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + mxc[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2);
				w[i] = interval(wdown[i], wup[i]); //% lower bound
			}
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<< 
//% The method is used in Krawczyk iteration to obtain the initital solution.
//% The method returns the midpoint inverse to decrease time complexity
//% of Krawczyk iteration.
//%----------------------------------------------------------------------------
HansenBliekRohnKRI(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w, //% interval vector solution
	aafrmatrix& V,
	aafrvector& v,
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), mxc(n), xast(n);
	dmatrix	Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = R * bc; //% midpoint solution
		mag(xc, mxc); //% abs of midpoint solution
		V = R * A;
		C = reduce(V); //% left pre-conditioning
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		v = R * b;
		z = reduce(v);
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (mxc + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - mxc[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + mxc[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2);
				w[i] = interval(wdown[i], wup[i]); //% lower bound
			}
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% The method combines the version with and without residual correction.
//%----------------------------------------------------------------------------
HansenBliekRohnComb(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xast(n), cmz(n), zr(n), zzr(n), wdown(n), wup(n), xcm(n), xc(n);
	dmatrix	R(n, n), Cr(n, n), M(n, n), I(n, n);
	ivector	z(n), zz(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		C = reduce(R*A);
		z = reduce(R*b);
		zz = reduce(R*(b - A * xc)); //% residual correction
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		rad(z, zr);
		rad(zz, zzr);
		mag(xc, xcm); //% xcm = |x^c|
		Unit(I);
		M = I - Cr;
		if (inv(M)) {
			xast = M * (xcm + zr);
			cmz = M * zzr;
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = M(i, i) * (xc[i] + xcm[i]) - xast[i];
				real v2 = ti * v1;
				wdown[i] = ISLIB_MAX(ISLIB_MIN(v1, v2), xc[i] - cmz[i]);
				v1 = M(i, i) * (xc[i] - xcm[i]) + xast[i];
				v2 = ti * v1;
				wup[i] = ISLIB_MIN(ISLIB_MAX(v1, v2), xc[i] + cmz[i]);
				w[i] = interval(wdown[i], wup[i]);
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnC <<<

int
//%----------------------------------------------------------------------------
//%  Ning & Kearfott method for solving PILS
//%  with affine matrix and affine right-hand vector
//%  Gives the same result as HBR method
//%----------------------------------------------------------------------------
NingKearfott(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	cmz(n), magz(n);
	dmatrix	Cm(n, n), S(n, n), I(n, n), R(n, n), Cr(n, n);
	ivector	z(n);
	imatrix	C(n, n), M(n, n);

	w = ivector(n);
	mid(A, R);
	if (inv(R)) {
		C = reduce(R * A); //% left pre-conditioning
		mig(C, Cm); //% Cm=<RA>
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		S = Cm;
		if (inv(S)) { //% Cm^{-1}
			z = reduce(R*b);
			mag(z, magz);
			dvector u = S * magz;
			real d, alpha, beta;
			for (int i = 0; i < n; i++) {
				d = S(i, i);
				alpha = Cm(i, i) - 1.0 / d;
				beta = u[i] / d - magz[i];
				w[i] = (z[i] + interval(-beta, beta)) / (C(i, i) + interval(-alpha, alpha));
			}
			return Success;
		}
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//%  Ning & Kearfott method for solving PILS
//%  with affine matrix and affine right-hand vector
//%  Gives the same result as HBR method
//%----------------------------------------------------------------------------
NingKearfottRC(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	x0(n), bc(n), cmz(n), magz(n);
	dmatrix	Cm(n, n), M(n, n), I(n, n), R(n, n), Cr(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	mid(A, R);
	if (inv(R)) {
		mid(b, bc);
		x0 = R * bc; //% midpoint solution
		C = reduce(R * A); //% left pre-conditioning
		mig(C, Cm); //% Cm=<RA>
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		M = Cm;
		if (inv(M)) { //% Cm^{-1}
			z = reduce(R*(b - A * x0));
			mag(z, magz);
			dvector u = M * magz;
			real d, alpha, beta;
			for (int i = 0; i < n; i++) {
				d = M(i, i);
				alpha = Cm(i, i) - 1.0 / d;
				beta = u[i] / d - magz[i];
				w[i] = x0[i] + (z[i] + interval(-beta, beta)) / (C(i, i) + interval(-alpha, alpha));
			}
			return Success;
		}
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//%  Ning & Kearfott method for solving PILS
//%  with affine matrix and affine right-hand vector
//%  Gives the same result as HBR method
//%----------------------------------------------------------------------------
NingKearfott(const imatrix& A, //% revised affine matrix
	const ivector& b, //% revised affine vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	cmz(n), magz(n);
	dmatrix	Cm(n, n), M(n, n), I(n, n), R(n, n), Cr(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	R = dmatrix(n, n);
	mid(A, R);
	if (inv(R)) {
		C = R * A; //% left pre-conditioning
		mig(C, Cm); //% Cm=<RA>
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		M = Cm;
		if (inv(M)) { //% Cm^{-1}
			z = R*b;
			mag(z, magz);
			dvector u = M * magz;
			real d, alpha, beta;
			for (int i = 0; i < n; i++) {
				d = M(i, i);
				alpha = Cm(i, i) - 1.0 / d;
				beta = u[i] / d - magz[i];
				w[i] = (z[i] + interval(-beta, beta)) / (C(i, i) + interval(-alpha, alpha));
			}
			return Success;
		}
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//%  Ning & Kearfott method for solving PILS
//%  with affine matrix and affine right-hand vector
//%  Gives the same result as HBR method
//%----------------------------------------------------------------------------
NingKearfottRC(const imatrix& A, //% revised affine matrix
	const ivector& b, //% revised affine vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	x0(n), bc(n), cmz(n), magz(n);
	dmatrix	Cm(n, n), M(n, n), I(n, n), R(n, n), Cr(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	R = dmatrix(n, n);
	mid(A, R);
	mid(b, bc);
	if (inv(R)) {
		x0 = R * bc; //% midpoint solution
		C = R * A; //% pleft pre-conditioning
		mig(C, Cm); //% Cm=<RA>
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		M = Cm;
		if (inv(M)) { //% Cm^{-1}
			z = R * (b - A * x0);
			mag(z, magz);
			dvector u = M * magz;
			real d, alpha, beta;
			for (int i = 0; i < n; i++) {
				d = M(i, i);
				alpha = Cm(i, i) - 1.0 / d;
				beta = u[i] / d - magz[i];
				w[i] = x0[i] + (z[i] + interval(-beta, beta)) / (C(i, i) + interval(-alpha, alpha));
			}
			return Success;
		}
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Improvement of the result obtained using Ning-Kerafott method
//%----------------------------------------------------------------------------
ErrorImpr(const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine vector
	ivector& w) //% interval solution vector
{
	//% This is probably not what they meant.
	//% Think it over!!! NGYUEN thesis strona 50
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), j = 0, err;
	dvector me(n);
	ivector e(n);
	aafrvector r(n), wf(n);

	err = NingKearfott(A, b, w);
	if (err != Success) return err;
	do {
		for (int i = 0; i < n; i++) {
			wf[i] = AAFR(w[i]);
		}
		err = NingKearfott(A, b - A * wf, e);
		if (err != Success) return err;
		mid(e, me);
		w = w + me;
		j++;
	} while (j < 20);
	return SingularMatrix; //% singular matrix
}

//%============================================================================
//%
//% Methods for solving parametric interval linear systems.
//% The methods use right preconditioning.
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% The method uses right preconditioning. The method returns the midpoint 
//% inverse. The method is used in Krawczyk method with right preconditioning.
//% The final result is premultiplied by R, because we need the solution,
//% which can be then improved by Krawczyk method using right preconditioning.
//% The final solution produced by Krawczyk method is premultiplied by R.
//%----------------------------------------------------------------------------
HansenBliekRohnII0(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w, //% interval vector solution
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n);
	dmatrix	Cr(n, n), I(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = bc; //% midpoint solution
		mag(xc, xcm); //% abs of midpoint solution
		C = reduce(A * R); //% right preconditioning
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		z = reduce(b);
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (xcm + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2);
				w[i] = interval(wdown[i], wup[i]); //% lower bound
			}
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnII0 <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% The method uses right preconditioning. The method returns the midpoint 
//% inverse.
//%----------------------------------------------------------------------------
HansenBliekRohnII(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n);
	ivector	z(n);
	dmatrix	R(n, n), Cr(n, n), I(n, n), M(n, n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint right-hand vector
		xc = bc; //% midpoint solution (we have A*R*y=b => midpoint solution is yc=bc)
		mag(xc, xcm); //% abs of midpoint solution
		C = reduce(A * R); //% right preconditioning
		rad(C, Cr); //% rad of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		z = reduce(b);
		rad(z, zr); //% radius of z
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			xast = M * (xcm + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2);
				w[i] = interval(wdown[i], wup[i]); //% lower bound
			}
			w = R * w;
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% >>> New operator by M. Hladik <<< 
//% for solving interval systems, with residual correction and preconditioning
//%----------------------------------------------------------------------------
NewOperatorMH(const aafrmatrix& A, //% interval matrix
	const aafrvector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	dvector	xc(n), magz(n), cmz(n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n), Br2(n, n);
	ivector	z(n), u(n);
	imatrix	B(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		B = reduce(R * A);
		z = reduce(R * b);
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		Br2 = Br * Br;
		if (rhoSpectral(Br) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			u = cmz * interval(-1.0, 1.0);
			real d, gamma, s;
			for (int i = 0; i < n; i++) {
				d = B(i, i).sup() / (1.0 - Br2(i, i));
				gamma = Bm(i, i) - 1.0 / d;
				s = 0.0;
				for (int j = 0; j < i; j++) {
					s += Br(i, j)*u[j].sup();
				}
				for (int j = i + i; j < n; j++) {
					s += Br(i, j)*u[j].sup();
				}
				w[i] = (z[i] + (s - gamma * u[i].inf()) * interval(-1.0, 1.0)) / (B(i, i) + gamma * interval(-1.0, 1.0));
			}
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> NewOperatorMH <<<

int
//%----------------------------------------------------------------------------
//% >>> New operator by M. Hladik <<< 
//% for solving interval systems, with residual correction and preconditioning
//%----------------------------------------------------------------------------
NewOperatorMHRC(const aafrmatrix& A, //% interval matrix
	const aafrvector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	interval ui(-1.0, 1.0);
	dvector	bc(n), xc(n), magz(n), cmz(n), x0(n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n), Br2(n, n);
	ivector	z(n), u(n);
	imatrix	B(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		mid(b, bc); //% midpoint vector
		x0 = R * bc;
		B = reduce(R * A);
		z = reduce(R * (b - A * x0));
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		Br2 = Br * Br;
		if (rhoSpectral(Br) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			u = cmz * ui;
			real d, gamma, s;
			for (int i = 0; i < n; i++) {
				d = B(i, i).sup() / (1.0 - Br2(i, i));
				gamma = Bm(i, i) - 1.0 / d;
				s = 0.0;
				for (int j = 0; j < i; j++) {
					s += Br(i, j)*u[j].sup();
				}
				for (int j = i + i; j < n; j++) {
					s += Br(i, j)*u[j].sup();
				}
				w[i] = x0[i] + (z[i] + (s - gamma * u[i].inf()) * ui) / (B(i, i) + gamma * ui);
			}
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> NewOperatorMHRC <<<

//%============================================================================
//%
//% Methods for solving parametric interval linear systems with multiple right-hand side
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn (HBR) method <<<
//% for solving systems with multiple right-hand side
//%----------------------------------------------------------------------------
HansenBliekRohnMRHS(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& b, //% revised affine right-hand matrix
	imatrix& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), m = w.num_cols();
	dmatrix	R(n, n), Cr(n, n), M(n, n), I(n, n);
	dmatrix	bc(n, m), Xast(n, m), Zr(n, m), y(n, m), wdown(n, m), wup(n, m), Xcm(n, m), Xc(n, m);
	imatrix	C(n, n), Z(n, m);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint vector
		Xc = R * bc; //% midpoint solution
		C = reduce(R * A); //% left pre-conditioning
		Z = reduce(R * b); //% left pre-conditioning
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		rad(Z, Zr);
		mag(Xc, Xcm); //% abs of midpoint solution
		Unit(I);
		M = I - Cr;
		if (inv(M)) {
			Xast = M * (Xcm + Zr);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					real ti = (1.0 / (2.0 * M(i, i) - 1.0));
					real v1 = M(i, i) * (Xc(i, j) + Xcm(i, j)) - Xast(i, j);
					real v2 = ti * v1;
					wdown(i, j) = ISLIB_MIN(v1, v2); //% lower bound
					v1 = M(i, i) * (Xc(i, j) - Xcm(i, j)) + Xast(i, j);
					v2 = ti * v1;
					wup(i, j) = ISLIB_MAX(v1, v2); //% upper bound
					w(i, j) = interval(wdown(i, j), wup(i, j));
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnMRHS <<<

//%============================================================================
//%
//% Methods for solving interval linear systems
//%
//%============================================================================

void
//%----------------------------------------------------------------------------
//% >>> BS versus HBR <<<
//% The method computes the difference between original BS and HBR bounds
//%----------------------------------------------------------------------------
BS_vs_HBR(const dmatrix& M,
	const dvector& xc,
	const dvector& xa,
	dvector& errl, //% difference between lower bounds
	dvector& erru) //% difference between upper bounds
{
	int	n = M.num_rows();
	real v1, v2;
	dvector xcm(n);

	mag(xc, xcm);
	errl = dvector(n);
	erru = dvector(n);
	for (int i = 0; i < n; i++) {
		v2 = 2.0 * (M(i, i) - 1.0) / (2.0 * M(i, i) - 1.0) * (xa[i] - xcm[i]);
		v1 = (M(i, i) - 1.0) * (xcm[i] + xc[i]);
		errl[i] = ISLIB_MIN(v1, v2); //% lower error
		v1 = (M(i, i) - 1.0) * (xcm[i] - xc[i]);
		erru[i] = ISLIB_MIN(v1, v2); //% upper error
	}
} //% >>> BS_vs_HBR <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn method <<<
//% for solving interval linear systems
//%----------------------------------------------------------------------------
HansenBliekRohn(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w, //% interval vector solution
	int err) //% 1 - print error
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n), errl, erru;
	dmatrix	R(n, n), Cr(n, n), M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		mag(xc, xcm); //% |x^c|
		z = R * b;
		C = R * A;
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		rad(z, zr); //% radius of z
		inf(C, M); //% equiv. to (I-rad([C]))
		if (inv(M)) { //% inverse of M
			xast = M * (xcm + zr); //% x*
			for (int i = 0; i < n; i++) {
				real ti = (1.0 / (2.0 * M(i, i) - 1.0));
				real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
				real v2 = ti * v1;
				wup[i] = ISLIB_MAX(v1, v2); //% upper bound
				v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
				v2 = ti * v1;
				wdown[i] = ISLIB_MIN(v1, v2);
				w[i] = interval(wdown[i], wup[i]); //% lower bound
			}
			if (err == 1) {
				BS_vs_HBR(M, xc, xast, errl, erru);
				for (int i = 0; i < n; i++) {
					std::cout << errl[i] << ",  " << erru[i] << std::endl;
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn method <<<
//% for solving interval systems. The method takes preconditioner R.
//% After preconditioning, the midpoint matrix is an identity matrix,
//% hence I-Cr = inf(C)
//%----------------------------------------------------------------------------
HansenBliekRohn(const dmatrix& R, //% pre-conditioner (midpoint inverse)
	const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval vector solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dvector	bc(n), xc(n), zr(n), wdown(n), wup(n), xcm(n), xast(n);
	dmatrix	M(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(b, bc); //% midpoint vector
	xc = R * bc; //% midpoint solution
	mag(xc, xcm); //% |x^c|
	z = R * b;
	C = R * A;
	rad(z, zr); //% zr = rad(z)
	inf(C, M); //% lower bound of C, inf(C) = I - rad(C)
	if (inv(M)) { //% M = (I-rad([C]))^{-1}
		xast = M * (xcm + zr); //% x^*
		for (int i = 0; i < n; i++) {
			real ti = (1.0 / (2.0 * M(i, i) - 1.0));
			real v1 = xast[i] + (xc[i] - xcm[i]) * M(i, i);
			real v2 = ti * v1;
			wup[i] = ISLIB_MAX(v1, v2); //% upper bound
			v1 = -xast[i] + (xc[i] + xcm[i]) * M(i, i);
			v2 = ti * v1;
			wdown[i] = ISLIB_MIN(v1, v2);
			w[i] = interval(wdown[i], wup[i]); //% lower bound
		}
		return Success;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohn <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn method <<< 
//% for solving interval systems, with residual correction
//%----------------------------------------------------------------------------
HansenBliekRohnRC(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	dvector	bc(n), xc(n);
	dmatrix	R(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		z = b - A * xc;
		err = HansenBliekRohn(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnRC <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn method <<< 
//% for solving interval systems, with preconditioning
//%----------------------------------------------------------------------------
HansenBliekRohnP(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dmatrix	R(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		C = R * A;
		z = R * b;
		return HansenBliekRohn(C, z, w, 0);
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnP <<<

int
//%----------------------------------------------------------------------------
//% >>> Hansen-Bliek-Rohn method <<< 
//% for solving interval systems, with residual correction and preconditioning
//%----------------------------------------------------------------------------
HansenBliekRohnRCP(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows();
	dmatrix	R(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		C = R * A;
		z = R * b;
		return HansenBliekRohnRC(C, z, w);
	}
	return SingularMatrix;
} //% >>> HansenBliekRohnRCP <<<

int
//%----------------------------------------------------------------------------
//% >>> New operator by M. Hladik <<< 
//% for solving interval systems, with residual correction and preconditioning
//%----------------------------------------------------------------------------
NewOperatorMH(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	interval ui(-1.0, 1.0);
	dvector	xc(n), magz(n), cmz(n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n), Br2(n, n);
	ivector	z(n), u(n);
	imatrix	B(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		B = R * A;
		z = R * b;
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		Br2 = Br * Br;
		if (rhoSpectral(Br) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			u = cmz * interval(-1.0, 1.0);
			real d, gamma, s;
			for (int i = 0; i < n; i++) {
				d = B(i, i).sup() / (1.0 - Br2(i, i));
				gamma = Bm(i, i) - 1.0 / d;
				s = 0.0;
				for (int j = 0; j < i; j++) {
					s += Br(i, j)*u[j].sup();
				}
				for (int j = i + i; j < n; j++) {
					s += Br(i, j)*u[j].sup();
				}
				w[i] = (z[i] + (s - gamma * u[i].inf()) * ui) / (B(i, i) + gamma * ui);
			}
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> New operator by M. Hladik <<< 

int
//%----------------------------------------------------------------------------
//% >>> New operator by M. Hladik <<< 
//% for solving interval systems, with residual correction and preconditioning
//%----------------------------------------------------------------------------
NewOperatorMHRC(const imatrix& A, //% interval matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	interval ui(-1.0, 1.0);
	dvector	bc(n), xc(n), magz(n), cmz(n), x0(n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n), Br2(n, n);
	ivector	z(n), u(n);
	imatrix	B(n, n);

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		mid(b, bc); //% midpoint vector
		x0 = R * bc;
		B = R * A;
		z = R * (b - A*x0);
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		Br2 = Br * Br;
		if (rhoSpectral(Br) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			u = cmz * ui;
			real d, gamma, s;
			for (int i = 0; i < n; i++) {
				d = B(i, i).sup() / (1.0 - Br2(i, i));
				gamma = Bm(i, i) - 1.0 / d;
				s = 0.0;
				for (int j = 0; j < i; j++) {
					s += Br(i, j)*u[j].sup();
				}
				for (int j = i + i; j < n; j++) {
					s += Br(i, j)*u[j].sup();
				}
				w[i] = x0[i] + (z[i] + (s - gamma * u[i].inf()) * ui) / (B(i, i) + gamma * ui);
			}
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> New operator by M. Hladik <<<