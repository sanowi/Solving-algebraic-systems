//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Direct methods for solving interval and parametric
//%     interval linear systems. The method requires that A is an H-matrix
//%     (strongly regular) and is based on the following formula:
//%     x^OI = xc + <RA>^{-1}|R(b-A*xc)|[-1,1].
//%     The method was described in the paper 
//%     "A method for outer interval solution...", Reliable Computing
//%
//%##########################################################################
//%

#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/randnum.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../truss/TrussStructure.h"
#include "../affine/aafr.h"

int
//%----------------------------------------------------------------------------
//% >>> Direct Method <<<
//% Method for solving parametric interval linear systems. The method uses
//% left pre-conditioning (to make the system matrix more tractable) and
//% residual correction of right-hand vector.
//%----------------------------------------------------------------------------
ParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand side vector
	ivector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	ivector	z(n);
	imatrix	B(n, n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n);
	dvector	bc(n), xc(n), magz(n), cmz(n);

	mid(A, R);  //% midpoint matrix
	mid(b, bc); //% midpoint vector
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		xc = R * bc; //% midpoint solution
		z = reduce(R*(b - A * xc)); //% left pre-conditioning and residual correction
		B = reduce(R*A); //% left pre-conditioning
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		if (rhoSpectral(Br) >= 1.0) { //% strong regularity verification
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			w = xc + cmz * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethod <<<

int 
//%----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems with multiple  
//% right-hand side.
//%----------------------------------------------------------------------------
ParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& B, //% revised right-hand side matrix
	imatrix& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	int m = B.num_cols();
	imatrix z(n, m), C(n, n);
	aafrmatrix Ax(n, m);
	dmatrix R(n, n), Cm(n, n), Cr(n, n);
	dmatrix bc(n, m), Xc(n, m), magz(n, m), cmz(n, m);

	w = imatrix(n, m);
	w.fill_in(0);
	mid(A, R); //% midpoint matrix
	mid(B, bc); //% midpoint right-hand matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		Xc = R * bc; //% midpoint solution
		z = reduce(R*(B - A * Xc)); //% left pre-conditioning and residual correction
		C = reduce(R*A); //% left pre-conditioning
		mig(C, Cm);
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) { //% strong regularity verification
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (inv(Cm)) {
			w = Xc + (Cm * magz) * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethod <<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems.
//% Additionally midpoint inverse, R, is returned.
//%----------------------------------------------------------------------------
ParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	ivector& w, //% p-solution (revised affine form)
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	ivector	z(n);
	imatrix	C(n, n);
	dmatrix	Cm(n, n), Cr(n, n);
	dvector	bc(n), x(n), magz(n), cmz(n);

	w.fill_in(0.0);
	mid(A, R); //% midpoint matrix
	mid(b, bc); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		x = R * bc; //% midpoint solution
		z = reduce(R * (b - A * x)); //% left pre-conditioning and residual correction
		C = reduce(R * A); //% left pre-conditioning
		mig(C, Cm);
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) { //% strong regularity verification
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = x + cmz * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethod<<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems with multiple 
//% right-handside. Additionally midpoint inverse, R, is returned.
//%-----------------------------------------------------------------------------
ParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& B, //% revised affine right-hand matrix
	imatrix& w, //% interval solution
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	int m = B.num_cols();
	imatrix z(n, m), C(n, n);
	aafrmatrix Ax(n, m);
	dmatrix Cm(n, n), Cr(n, n), S(n, n);
	dmatrix bc(n, m), Xc(n, m), magz(n, m), cmz(n, m);

	w = imatrix(n, m);
	w.fill_in(0);
	R = dmatrix(n, n);
	mid(A, R); //% midpoint matrix
	mid(B, bc); //% midpoint right-hand matrix
	if (inv(R)) { //% midpoint inverse
		Xc = R * bc; // midpoint solution
		z = reduce(R*(B - A * Xc)); //% left pre-conditioning and residual correction
		C = reduce(R * A); //% left pre-conditioning
		mig(C, Cm); //% Cm=<C>
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) { //% strong regularity verification
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (inv(Cm)) {
			Xc + (Cm * magz) * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethod<<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems. The method uses
//% right preconditioning and residual correction. Additionally midpoint 
//% inverse is returned.
//%-----------------------------------------------------------------------------
ParametricDirectMethodPII(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	ivector& w, //% p-solution (revised affine form)
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	ivector z(n);
	imatrix C(n, n);
	dmatrix Cm(n, n), Cr(n, n);
	dvector bc(n), yc(n), magz(n), cmz(n);
	aafrmatrix B(n, n);

	w.fill_in(0.0);
	mid(A, R); //% revised affine matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% revised right-hand vector
		yc = bc; //% midpoint solution of the systen By=b (mid(B)=I => xc=bc)
		B = A * R; //% right pre-conditioning
		z = reduce(b - B * yc); //% residual correction
		C = reduce(B);
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = yc + cmz * interval(-1.0, 1.0); //% solution to the system By=b
			return Success;
		}
	}
	return SingularMatrix; //% singular midpoint matrix => A singular
} //% >>> ParametricDirectMethodII <<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<< 
//% Method for solving parametric interval linear systems. The method uses
//% right preconditioning and residual correction. Additionally midpoint
//% inverse is returned.
//%-----------------------------------------------------------------------------
ParametricDirectMethodP(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	ivector& w, //% p-solution (revised affine form)
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	ivector z(n);
	imatrix C(n, n);
	dmatrix Cm(n, n), Cr(n, n);
	dvector bc(n), yc(n), magz(n), cmz(n);
	aafrmatrix B(n, n);

	w.fill_in(0.0);
	mid(A, R); //% midpoint matrix
	mid(b, bc); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		yc = bc; //% first we pre-condition with midpoint inverse => mid(B)=I and xc = bc
		B = A * R; //% right pre-conditioning
		z = reduce(b - B * yc); //% residual correction
		C = reduce(B);
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) { 
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = R * (yc + cmz * interval(-1.0, 1.0));
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethodP <<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems. The method uses
//% right preconditioning and residual correction. Additionally midpoint 
//% inverse is returned.
//%-----------------------------------------------------------------------------
ParametricDirectMethodPIII(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	const dmatrix& L, //% left preonditioner
	const dmatrix& R, //% right preconditioner
	ivector& w) //% p-solution (revised affine form)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector bb(n);
	ivector z(n);
	imatrix C(n, n);
	dmatrix Ac(n, n), Cm(n, n), Cr(n, n);
	dvector bc(n), xc(n), magz(n), cmz(n);
	aafrmatrix B(n, n);

	w.fill_in(0.0);
	mid(A, Ac);
	mid(b, bc); //% revised right-hand vector
	if (inv(Ac)) { //% midpoint inverse
		B = L * A * R; //% left and right pre-conditioning
		xc = Ac * bc; //% midpoint solution  - can be obtained using e.g. Gauss elimination
		z = reduce(L*(b - A * xc));
		C = reduce(B);
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = xc + cmz * interval(-1.0, 1.0);
			return Success;
		}
	}
} //% >>> ParametricDirectMethodIII <<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems. The method uses
//% right preconditioning and residual correction. Additionally midpoint 
//% inverse is returned.
//%-----------------------------------------------------------------------------
ParametricDirectMethodPIII2(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	const dmatrix& L, //% left preonditioner
	const dmatrix& R, //% right preconditioner
	ivector& w) //% p-solution (revised affine form)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector bb(n);
	ivector z(n);
	imatrix C(n, n);
	dmatrix Bc(n, n), Cm(n, n), Cr(n, n);
	dvector bc(n), yc(n), magz(n), cmz(n);
	aafrmatrix B(n, n), D(n, n);

	w.fill_in(0.0);
	B = A * R;
	mid(B, Bc);
	mid(b, bc); //% revised right-hand vector
	if (inv(Bc)) { //% midpoint inverse
		D = L * B; //% left and right pre-conditioning
		yc = Bc * bc; //% midpoint solution  - can be obtained using e.g. Gauss elimination
		z = reduce(L*(b - B * yc));
		C = reduce(D);
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = yc + cmz * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethodIII2 <<<

int
//%-----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving parametric interval linear systems. The method uses
//% left and right preconditioning and residual correction. 
//%-----------------------------------------------------------------------------
ParametricDirectMethodPP(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised right-hand vector
	const dmatrix& L, //% left preonditioner
	const dmatrix& R, //% right preconditioner
	ivector& w) //% p-solution (revised affine form)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrvector bb(n);
	ivector z(n);
	imatrix C(n, n);
	dmatrix Ac(n, n), Cm(n, n), Cr(n, n);
	dvector bc(n), xc(n), magz(n), cmz(n);
	aafrmatrix B(n, n);

	w.fill_in(0.0);
	mid(A, Ac);
	mid(b, bc); //% revised right-hand vector
	if (inv(Ac)) { //% midpoint inverse
		xc = Ac * bc; //% midpoint solution
		B = L * A * R; //% left and right pre-conditioning
		z = reduce(L*(b - A * xc));
		C = reduce(B);
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		rad(C, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			w = xc + R * (cmz * interval(-1.0, 1.0));
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethodII <<<

int
//%----------------------------------------------------------------------------
//% >>> Parametric Direct Method <<<
//% Method for solving interval linear systems. The method uses
//% left pre-conditioning in order to make the system matrix more tractable.
//% The final solution is given as an interval vector.
//%----------------------------------------------------------------------------
ParametricDirectMethodM(const imatrix& A, //% revised affine system matrix
	const ivector& b, //% revised right-hand side vector
	ivector& w) //% interval solution
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	ivector	z(n);
	imatrix	B(n, n);
	dmatrix	R(n, n), Bm(n, n), Br(n, n);
	dvector	bm(n), xc(n), magz(n), cmz(n);

	mid(A, R);  //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		xc = R * bm; //% midpoint solution
		B = R * A;
		z = R * (b - A * xc);
		mig(B, Bm); //% Bm=<B> (Ostrowski comparison matrix)
		rad(B, Br);
		if (rhoSpectral(Br) >= 1.0) {
			return NotStronglyRegular;
		}
		mag(z, magz); //% magz=|z|
		if (gw_gausse(Bm, magz, cmz)) {
			w = xc + cmz * interval(-1.0, 1.0);
			return Success;
		}
	}
	return SingularMatrix;
} //% >>> ParametricDirectMethodP <<<

int
//%-----------------------------------------------------------------------------
//% >>> Kolev's Parametric Direct method <<< 
//% Modified Parametric Direct method, which prodeuces p-solution
//% (described by Kolev)
//%-----------------------------------------------------------------------------
KParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& w) //% p-solution (revised affine vector)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix Cinv(n, n);
	aafrvector z(n), zf(n);
	dmatrix R(n, n), Cr(n, n), I(n, n), M(n, n);
	dvector bm(n), xc(n);
	imatrix Ci(n, n), H(n, n);
	cvector idx(1);
	dvector cfs(1);

	idx[0] = 0;
	cfs[0] = 0.0;
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bm); //% midpoint vector
		xc = R * bm; //% midpoint solution
		z = R * (b - A * xc); //% left pre-conditioning and residual correction
		Ci = reduce(R*A); //% left pre-conditioning
		rad(Ci, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M computed using Rohn's formula
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i == j) {
						H(i, j) = interval(M(i, j) / (2.0 * M(i, j) - 1.0), M(i, j));
					}
					else {
						H(i, j) = interval(-M(i, j), M(i, j));
					}
				}
			}
			//% convert H into affine matrix
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					cfs[0] = H(i, j).mid();
					Cinv(i, j) = AAFR(idx, cfs, H(i, j).rad());
				}
			}
			zf = Cinv * z;
			for (int i = 0; i < n; i++) {
				w[i] = xc[i] + zf[i];
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> KParametricDirectMethod <<<

int
//%------------------------------------------------------------------------------
//% >>> Kolev's Parametric Direct method <<< 
//% Parametric version of Direct method (ISM).
//% The method uses right preconditioning
//%-----------------------------------------------------------------------------
KParametricDirectMethodII(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& w) //% p-solution (revised affine vector)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix B(n, n), Cinv(n, n);
	aafrvector z(n), zf(n);
	dmatrix R(n, n), Cr(n, n), I(n, n), M(n, n);
	dvector bc(n), xc(n);
	imatrix Ci(n, n), H(n, n);
	cvector idx(1);
	dvector cfs(1);

	idx[0] = 0;
	cfs[0] = 0.0;
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint vector
		B = A * R; //% right pre-conditioning
		xc = R * bc; //% midpoint solution
		z = b - A * xc; //% residual correction
		Ci = reduce(B);
		rad(Ci, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		Unit(I);
		M = I - Cr;
		if (inv(M)) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i == j) {
						H(i, j) = interval(M(i, j) / (2.0 * M(i, j) - 1.0), M(i, j));
					}
					else {
						H(i, j) = interval(-M(i, j), M(i, j));
					}
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					cfs[0] = H(i, j).mid();
					Cinv(i, j) = AAFR(idx, cfs, H(i, j).rad());
				}
			}
			zf = R * Cinv * z;
			for (int i = 0; i < n; i++) {
				w[i] = xc[i] + zf[i];
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> KParametricDirectMethodII <<<

int
//%-----------------------------------------------------------------------------
//% >>> Kolev's Parametric Direct method <<< 
//% Parametric version of Direct method (ISM).
//% Method solves systems with multiple right-hand side.
//%-----------------------------------------------------------------------------
KParametricDirectMethod(const aafrmatrix& A, //% revised affine system matrix
	const aafrmatrix& b, //% revised affine right-hand side matrix
	aafrmatrix& w) //% p-solution (revised affine matrix)
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = b.num_cols();
	aafrmatrix AA(n, n), Cinv(n, n);
	aafrmatrix bb(n, m), z(n, m), zf(n, m);
	dmatrix R(n, n), Cr(n, n), I(n, n), M(n, n);
	dmatrix bm(n, m), Xc(n, m);
	imatrix Ci(n, n), H(n, n);
	cvector idx(1);
	dvector cfs(1);
	AAFR s;

	idx[0] = 0;
	cfs[0] = 0.0;

	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bm); //% midpoint vector
		Xc = R * bm; //% midpoint solution
		z = R * (b - A * Xc); //% left pre-conditioning and residual correction
		Ci = reduce(R * A); //% left pre-conditioning
		rad(Ci, Cr);
		if (rhoSpectral(Cr) >= 1.0) {
			return NotStronglyRegular;
		}
		Unit(I);
		M = I - Cr;
		if (inv(M)) {
			//% H is the (Rohn) inverse of  interval matrix
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i == j) {
						H(i, j) = interval(M(i, j) / (2.0 * M(i, j) - 1.0), M(i, j));
					}
					else {
						H(i, j) = interval(-M(i, j), M(i, j));
					}
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					cfs[0] = H(i, j).mid();
					Cinv(i, j) = AAFR(idx, cfs, H(i, j).rad());
				}
			}
			zf = Cinv * z;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					w(i, j) = Xc(i, j) + zf(i, j);
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
} //% >>> KParametricDirectMethod <<<

void
//%-----------------------------------------------------------------------------
//% >>> OUSystem <<<
//% Creates augmented system from over- or under-determined systems. The 
//% size of the new system is n + m.
//%-----------------------------------------------------------------------------
OUSystem(const aafrmatrix& A, //% revised affine system matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrmatrix& AA, //% augmented revised affine matrix 
	aafrvector& bb) //% augmented revised affine right-hand vector
{
	int	n = A.num_rows(), m = A.num_cols(), nm = n + m;
	dvector	x(nm), bm(nm), cmz(nm), magz(nm);
	dmatrix	R(nm, nm), Cm(nm, nm);
	ivector	z(nm);
	imatrix	C(nm, nm);

	AA = aafrmatrix(nm, nm);
	bb = aafrvector(nm);
	nm = n + m;
	AA.fill_in(0.0);
	bb.fill_in(0.0);
	if (n > m) { //% overdetermined
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AA(i, j) = A(i, j);
			}
			AA(i, i + m) = -1.0;
			bb[i] = b[i];
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				AA(i + n, j + m) = A(j, i);
			}
		}
	}
	else { //% underdetermined
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				AA(i, j) = A(j, i);
			}
			AA(i, i + n) = -1.0;
			bb[i + m] = b[i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AA(i + m, j + n) = A(i, j);
			}
		}
	}
} //% >>> OUSystem <<<