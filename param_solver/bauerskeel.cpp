//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Bauer-Skeel method
//%
//%##########################################################################
//%

#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/direct_solvers.h"
#include "../param_solver/old_methods.h"
#include "../truss/TrussStructure.h"
#include "../examples/pexamples.h"
#include "../param_solver/bauerskeel.h"

//%============================================================================
//%
//% Methods for solving parametric interval linear systems
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% * without residual correction
//% * preconditioned with the midpoint inverse
//% * after preconditoning the formulas from the Rohn's paper are used.
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand side vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector bm(n), xc(n), zr(n), xast(n), xcm(n);
	dmatrix R(n, n), M(n, n), I(n, n), Cr(n, n);
	ivector z(n);
	imatrix C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse - preconditioner
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% tilde x - midpoint solution
		C = reduce(R*A);
		rad(C, Cr); //% radius of C
		rs = rhoSpectral(Cr);
		if (rs >= 1) {
			return NotHMatrix; //% the necessary condition is not fulfilled
		}
		z = reduce(R*b);
		rad(z, zr); //% radius of z
		mag(xc, xcm); //% abs of the midpoint solution
		Unit(I);
		M = I - Cr;
		if (gw_gausse(M, Cr*xcm + zr, xast)) { //% x*
			w = xc + xast * interval(-1.0, 1.0);
			return Success;
		}
		return SingularMatrix;
	}
	return SingularMatrix; //% singular matrix
} //% >>> BauerSkeel <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% * with residual correction
//% * preconditioned with the midpoint inverse
//% * after preconditoning the formulas from the Rohn's paper are used.
//%----------------------------------------------------------------------------
BauerSkeelRC(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector z(n), zdm(n);
	aafrmatrix C(n, n);	
	dvector bm(n), xc(n);
	dmatrix R(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		xc = R * bm; //% midpoint solution
		z = b - A * xc; //% residual correction
		err = BauerSkeel(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix; //% singular matrix
} //% >>> BauerSkeelRC <<<

int
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction and with additional preconditioning. Here left
//% pre-conditioning with midpoint inverse (which is optimal) is used, however, 
//% any other nonsingular matrix R can be used.
//%----------------------------------------------------------------------------
BauerSkeelP(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% refived affine right-hand side vector
	int r, //% residual correction
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return  NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector z(n);
	aafrmatrix C(n, n);	
	dvector bm(n), xc(n);
	dmatrix R(n, n);

	w = ivector(n);
	mid(A, R);
	if (inv(R)) { //% midpoint inverse - preconditioner
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% midpoint solution
		C = R * A; //% left pre-conditioning
		if (r == 0) {
			z = R * b; //% left pre-conditioning
			return BauerSkeel(C, z, w);
		}
		z = R * (b - A * xc);
		err = BauerSkeel(C, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix; //% singular matrix
} //% >>> BauerSkeelP method <<<

int
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction. Implementation corresponds to left
//% pre-conditioning with the midpoint inverse. The method uses precomputed 
//% midpoint inverse.
//%----------------------------------------------------------------------------
BauerSkeel(const dmatrix& R, //% pre-conditioner (midpoint inverse)
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% refived affine right-hand side vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector bm(n), xc(n), zr(n), xast(n), xcm(n);
	dmatrix M(n, n), I(n, n), Cr(n, n);
	ivector z(n);
	imatrix C(n, n);

	w = ivector(n);
	mid(b, bm); //% midpoint right-hand vector
	xc = R * bm; //% midpoint solution
	C = reduce(R*A);
	rad(C, Cr); //% Cr = rad(C)
	rs = rhoSpectral(Cr);
	if (rs >= 1) {
		return NotHMatrix; //% the necessary condition is not fulfilled
	}
	z = reduce(R*b);
	rad(z, zr); //% zr = rad(z)
	mag(xc, xcm); //% abs of the midpoint solution
	Unit(I);
	M = I - Cr;
	if (gw_gausse(M, Cr*xcm + zr, xast)) { //% x*
		w = xc + xast * interval(-1.0, 1.0);
		return Success;
	}
	return SingularMatrix;
} //% >>> BauerSkeel <<<

int
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction. The formulas correspond to he original 
//% version of the BS method described by Rohn.
//%----------------------------------------------------------------------------
BauerSkeel1(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% refived affine right-hand side vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector	bm(n), br(n), xc(n), zr(n), xast(n), xcm(n), bb(n);
	dmatrix	R(n, n), M(n, n), I(n, n), Cr(n, n), Rm(n, n), Ar(n, n);
	
	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse - preconditioner
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% midpoint solution
		mag(R, Rm); //% abs of midpoint inverse
		rad(reduce(b), br); //% radius of b
		rad(reduce(A), Ar); //% radius of A
		Cr = Rm * Ar;
		rs = rhoSpectral(Cr);
		if (rs >= 1) {
			return NotHMatrix; //% the necessary condition is not fulfilled
		}
		mag(xc, xcm); //% abs of midpoint solution
		bb = Rm * (Ar * xcm + br);
		Unit(I);
		M = I - Cr;
		if (gw_gausse(M, bb, xast)) { //% x^*
			w = xc + xast * interval(-1.0, 1.0);
			return Success;
		}
		return SingularMatrix;
	}
	return SingularMatrix; //% singular matrix
} //% >>> BauerSkeel1 <<<

int
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction 
//% corresponds to the method which is based on the original version of BS.
//%----------------------------------------------------------------------------
BauerSkeel1(const dmatrix& R, //% pre-conditioner (midpoint inverse)
	const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% refived affine right-hand side vector
	ivector& w) //% interval solution vector

{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector	bm(n), br(n), xc(n), zr(n), xast(n), xcm(n), bb(n);
	dmatrix	M(n, n), I(n, n), Cr(n, n), Rm(n, n), Ar(n, n);

	w = ivector(n);
	mid(b, bm); //% midpoint right-hand vector
	xc = R * bm; //% midpoint solution
	mag(R, Rm); //% abs of midpoint inverse
	rad(reduce(b), br); //% radius of b
	rad(reduce(A), Ar); //% radius of A
	mag(xc, xcm); //% abs of midpoint solution
	bb = Rm * (Ar * xcm + br);
	Cr = Rm * Ar;
	rs = rhoSpectral(Cr);
	if (rs >= 1) {
		return NotHMatrix; //% the necessary condition is not fulfilled
	}
	Unit(I);
	M = I - Cr;
	if (gw_gausse(M, bb, xast)) { //% x*
		w = xc + xast * interval(-1.0, 1.0);
		return Success;
	}
	return SingularMatrix;
} //% >>> BauerSkeel1 <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with residual correction
//% corresponds to the method which is based on the original version of BS.
//%----------------------------------------------------------------------------
BauerSkeel1RC(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), err;
	aafrvector z(n), zdm(n);
	aafrmatrix C(n, n);
	dvector bm(n), xc(n);
	dmatrix R(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		xc = R * bm; //% midpoint solution
		z = b - A * xc; //% residual correction
		err = BauerSkeel1(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix;
} //% >>> BauserSkeel1RC <<<

//%============================================================================
//%
//% Methods for solving parametric interval linear systems with multiple right-hand side
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning nd residual correction for solving systems with 
//% multiple right-hand side. The method corresponds to BauerSkeelRC method.
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& B, //% revised affine right-hand matrix
	imatrix& W) //% revised affine matrix solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = W.num_cols();
	real rs;
	dmatrix	R(n, n), M(n, n), I(n, n), Cr(n, n), Bm(n, m);
	dmatrix	Xc(n, m), Zr(n, m), Xcm(n, m);
	imatrix	Z(n, n), C(n, n);

	W = imatrix(n, m);
	mid(A, R); //% midpoint matrix
	mid(B, Bm); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		Xc = R * Bm; //% midpoint solution
		C = reduce(R * A); //% left pre-conditioning
		rad(C, Cr); //% radius of C
		rs = rhoSpectral(Cr);
		if (rs >= 1) {
			return NotHMatrix; //% the necessary condition is not fulfilled
		}
		Z = reduce(R * (B - A * Xc)); //% left pre-conditioning and residual correction
		rad(Z, Zr); //% radius of Z
		mag(Xc, Xcm); //% abs of the midpoint solution
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			W = Xc + (M * Zr) * interval(-1.0, 1.0);
			return Success;
		}
		return SingularMatrix;
	}
	return SingularMatrix;
} //% >>> BauserSkeel <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning and residual correction
//% for solving systems with multiple right-hand side. Additionally,
//% the midpoint inverse is returned
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& B, //% revised affine right-hand matrix
	imatrix& W, //% revised affine matrix solution
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = W.num_cols();
	real rs;
	dmatrix	M(n, n), I(n, n), Cr(n, n);
	dmatrix	Bm(n, m), Xc(n, m), Zr(n, m), Xcm(n, m), Xast(n, m);
	imatrix	C(n, n), Z(n, m);

	W = imatrix(n, m);
	R = dmatrix(n, n);
	mid(A, R); //% midpoint matrix
	mid(B, Bm); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		Xc = R * Bm; //% midpoint solution
		C = reduce(R * A); //% left pre-conditioning
		Z = reduce(R * (B - A * Xc)); //% left pre-conditioning and residual correction
		rad(Z, Zr); //% radius of Z
		rad(C, Cr); //% radius of C
		rs = rhoSpectral(Cr);
		if (rs >= 1) {
			return NotHMatrix; //% the necessary condition is not fulfilled
		}
		mag(Xc, Xcm); //% abs of the midpoint solution
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			Xast = M * Zr;
			W = Xc + Xast * interval(-1.0, 1.0);
			return Success;
		}
		return SingularMatrix;
	}
	return SingularMatrix;
} //% >>> BauserSkeel <<<

//%============================================================================
//%
//% Methods for solving interval linear systems
//%
//%============================================================================

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<< 
//% for solving interval linear systems standard version (after Rohn)
//%----------------------------------------------------------------------------
BauerSkeel(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector	bm(n), xc(n), zr(n), xcm(n);
	dmatrix	R(n, n), M(n, n), Cr(n, n), Rm(n, n), Ar(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint right-hand vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		xc = R * bm; //% midpoint solution
		C = R * A; //% left pre-conditioning
		inf(C, M); //% M=inf([A])=I-|R|*rad([A])
		rad(C, Cr); //% Cr=rad([C])=|R|*rad([A])
		rs = rhoSpectral(Cr);
		if (rs >= 1) {
			return NotHMatrix; //% the necessary condition is not fulfilled
		}
		z = R * b; //% left pre-conditioning
		rad(z, zr); //% zr=rad([z])=|R|*rad([b])
		mag(xc, xcm); //% abs of midpoint solution
		if (inv(M)) { //% inverse of M
			w = xc + M * (Cr * xcm + zr) * interval(-1.0, 1.0);
			return Success;
		}
		return SingularMatrix;
	}
	return SingularMatrix;
} //% >>> BauerSkeel <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning for solving interval linear systems. The midpoint
//% inverse is used here, however any other preconditioner can be used.
//%----------------------------------------------------------------------------
BauerSkeelP(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	dvector bm(n), xc(n), zr(n), xcm(n);
	dmatrix R(n, n), M(n, n), Cr(n, n);
	ivector z(n);
	imatrix C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = R * A; //% left pre-conditioning
		z = R * b; //% left pre-conditioning
		return BauerSkeel(C, z, w);
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<< 
//% for solving interval linear systems standard version (after Rohn).  
//% The method takes precalculated preconditioner.
//%----------------------------------------------------------------------------
BauerSkeel(const dmatrix& R, //% pre-conditioner
	const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	real rs;
	dvector	bm(n), xc(n), zr(n), xcm(n);
	dmatrix	M(n, n), Cr(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	w = ivector(n);
	mid(b, bm); //% midpoint right-hand vector
	xc = R * bm; //% midpoint solution
	C = R * A; //% left pre-conditioning
	z = R * b; //% left pre-conditioning
	inf(C, M); //% M=inf([A])=I-|R|*rad([A])
	rad(C, Cr); //% radius of C
	rs = rhoSpectral(Cr);
	if (rs >= 1) {
		return NotHMatrix; //% the necessary condition is not fulfilled
	}
	rad(z, zr); //% radius of z
	mag(xc, xcm); //% abs of midpoint solution
	if (inv(M)) { //% inverse of M
		w = xc + M * (Cr*xcm + zr) * interval(-1.0, 1.0);
		return Success;
	}
	return SingularMatrix;
} //% >>> BauerSkeel <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with residual correction
//% for solving interval linear systems
//%----------------------------------------------------------------------------
BauerSkeelRC(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int	n = A.num_rows(), err;
	dvector bm(n), xc(n);
	dmatrix R(n, n);
	ivector z(n);
	imatrix C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse - preconditioner
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% midpoint solution
		z = b - A * xc; //% residual correction
		err = BauerSkeel(R, A, z, w);
		if (err == Success) {
			w = xc + w;
		}
		return err;
	}
	return SingularMatrix;
} //% >>> BauerSkeelRC <<<

int
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning for solving interval linear systems. The midpoint
//% inverse is used here, however any other preconditioner can be used.
//%----------------------------------------------------------------------------
BauerSkeelPRC(const imatrix& A, //% interval system matrix
	const ivector& b, //% interval right-hand vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	dvector bm(n), xc(n), zr(n), xcm(n), zc(n);
	dmatrix R(n, n), M(n, n), Cr(n, n);
	ivector		z(n), zz(n);
	imatrix		C(n, n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = R * A; //% left pre-conditioning
		z = R * b; //% left pre-conditioning
		return BauerSkeelRC(C, z, w);
	}
	return SingularMatrix;
} //% >>> BauerSkeelPRC <<<