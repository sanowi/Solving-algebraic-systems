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
#include "../examples/acmExamples.h"
#include "../param_solver/bauerskeel.h"

//%============================================================================
//%
//% Methods for solving paramteric interval linear systems
//%
//%============================================================================

bool
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% The system is first preconditioned with the midpoint inverse and then 
//% the formulas from the Rohn's paper are used. If r != 0, then residual
//% correction is applied.
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector&			b, //% refived affine right-hand side vector
	int							r, //% residual correction (0 - no, 1 - yes)
	ivector&					w) //% interval solution
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), M(n, n), I(n, n), Cr(n, n);
	dvector		bm(n), xc(n), zr(n), xast(n), xcm(n), bb(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% tilde x - midpoint solution
		C = reduce(R*A); //% preconditioning
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1) {
			return false; //% the necessary condition is not fulfilled
		}
		mag(xc, xcm); //% abs of midpoint solution
		if (r == 0) {
			z = reduce(R*b); //% preconditioning
			rad(z, zr); //% radius of z
			bb = Cr * xcm + zr;
		}
		else {
			z = reduce(R*(b - A * xc));
			rad(z, zr); //% radius of z
			bb = zr;
		}
		Unit(I);
		M = I - Cr;
		if (gw_gausse(M, bb, xast)) { //% x*
			w = xc + xast * interval(-1.0, 1.0);
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauerSkeel <<<

bool
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction and with additional preconditioning. Here
//% preconditioning with midpoint inverse (which is optimal) is used, however, 
//% any other nonsingular matrix R can be used.
//%----------------------------------------------------------------------------
BauerSkeelP(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector&			b, //% refived affine right-hand side vector
	int							r, //% residual correction
	ivector&					w) //% interval vector solution
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	aafrvector	z(n);
	aafrmatrix	C(n, n);
	dmatrix		R(n, n);
	dvector		bm(n), xc(n);

	w = ivector(n);
	mid(A, R);
	mid(b, bm);
	if (inv(R)) {
		xc = R * bm;
		C = R * A; //% preconditioning
		if (r == 0) {
			z = R * b; //% preconditioning
			return BauerSkeel(C, z, 0, w);
		}
		z = R * (b - A * xc);
		if (BauerSkeel(C, z, 0, w)) {
			w = xc + w;
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauerSkeelP method <<<

bool
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction. The formulas correspond to he original 
//% version of the BS method described by Rohn.
//%----------------------------------------------------------------------------
BauerSkeel1(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector&			b, //% refived affine right-hand side vector
	ivector&					w) //% interval solution
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), M(n, n), I(n, n), Cr(n, n), Rm(n, n), Ar(n, n);
	dvector		bm(n), br(n), xc(n), zr(n), xast(n), xcm(n), bb(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bm); //% midpoint right-hand vector
		xc = R * bm; //% midpoint solution
		mag(R, Rm); //% abs of midpoint inverse
		rad(reduce(b), br); //% radius of b
		rad(reduce(A), Ar); //% radius of A
		Cr = Rm * Ar;
		if (rhoSpectral(Cr) >= 1) {
			return false; //% the necessary condition is not fulfilled
		}
		mag(xc, xcm); //% abs of midpoint solution
		bb = Rm * (Ar * xcm + br);
		Unit(I);
		M = I - Cr;
		if (gw_gausse(M, bb, xast)) { //% x^*
			w = xc + xast * interval(-1.0, 1.0);
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauerSkeel1 <<<

bool
//%----------------------------------------------------------------------------
//% >>> Bauer-Skeel method <<<
//% without residual correction, the formulas correspond to the original 
//% version of the BS method described by Rohn. The method uses precomputed
//% midpoint inverse.
//%----------------------------------------------------------------------------
BauerSkeel1(const dmatrix&	R, //% precomputed midpoint inverse
	const aafrmatrix&		A, //% revised affine system matrix
	const aafrvector&		b, //% refived affine right-hand side vector
	ivector&				w) //% interval solution

{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		M(n, n), I(n, n), Cr(n, n), Rm(n, n), Ar(n, n);
	dvector		bm(n), br(n), xc(n), zr(n), xast(n), xcm(n), bb(n);

	w = ivector(n);
	mid(b, bm); //% midpoint right-hand vector
	xc = R * bm; //% midpoint solution
	mag(R, Rm); //% abs of midpoint inverse
	rad(reduce(b), br); //% radius of b
	rad(reduce(A), Ar); //% radius of A
	mag(xc, xcm); //% abs of midpoint solution
	bb = Rm * (Ar * xcm + br);
	Cr = Rm * Ar;
	if (rhoSpectral(Cr) >= 1) {
		return false; //% the necessary condition is not fulfilled
	}
	Unit(I);
	M = I - Cr;
	if (gw_gausse(M, bb, xast)) { //% x^*
		w = xc + xast * interval(-1.0, 1.0);
		return true;
	}
	return false; //% singular matrix
} //% >>> BauerSkeel1 <<<

bool
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with residual correction
//% corresponds to the method which is based on the original version of BS.
//%----------------------------------------------------------------------------
BauerSkeel1RC(const aafrmatrix&	A, //% revised affine system matrix
	const aafrvector&			b, //% revised affine right-hand vector
	ivector&					w) //% interval vector solution
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	aafrvector	z(n), zdm(n);
	aafrmatrix	C(n, n);
	dmatrix		R(n, n);
	dvector		bm(n), xc(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		xc = R * bm; //% midpoint solution
		z = b - A * xc; //% residual correction
		if (BauerSkeel1(R, A, z, w)) {
			w = xc + w;
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauserSkeel1RC <<<

//%============================================================================
//%
//% Methods for solving interval linear systems with multiple right-hand side
//%
//%============================================================================

bool
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning nd residual correction for solving systems with 
//% multiple right-hand side. The method corresponds to BauerSkeelRC method.
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix&	A, //% revised affine matrix
	const aafrmatrix&			B, //% revised affine right-hand matrix
	imatrix&					W) //% revised affine matrix solution
{
	if (!A.is_square()) return false;
	int			n = A.num_rows(), m = W.num_cols();
	imatrix		Z(n, n);
	imatrix		C(n, n);
	dmatrix		R(n, n), M(n, n), I(n, n), Cr(n, n), Bm(n, m);
	dmatrix		Xc(n, m), Zr(n, m), Xcm(n, m);

	W = imatrix(n, m);
	mid(A, R); //% midpoint matrix
	mid(B, Bm); //% midpoint vector
	if (inv(R)) {
		Xc = R * Bm; //% midpoint solution
		C = reduce(R * A); //% preconditioning
		Z = reduce(R * (B - A * Xc)); //% preconditioning and residual correction
		rad(Z, Zr); //% radius of Z
		rad(C, Cr); //% radius of C
		if (rhoSpectral(Cr) >= 1) {
			return false; //% the necessary condition is not fulfilled
		}
		mag(Xc, Xcm); //% abs of the midpoint solution
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			W = Xc + (M * Zr) * interval(-1.0, 1.0);
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauserSkeel <<<

bool
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning and residual correction
//% for solving systems with multiple right-hand side. Additionally,
//% the midpoint inverse is returned
//%----------------------------------------------------------------------------
BauerSkeel(const aafrmatrix&	A, //% revised affine matrix
	const aafrmatrix&			B, //% revised affine right-hand matrix
	imatrix&					W, //% revised affine matrix solution
	dmatrix&					R) //% midpoint inverse
{
	if (!A.is_square()) return false;
	int			n = A.num_rows(), m = W.num_cols();
	imatrix		C(n, n), Z(n, m);
	dmatrix		M(n, n), I(n, n), Cr(n, n);
	dmatrix		Bm(n, m), Xc(n, m), Zr(n, m), Xcm(n, m), Xast(n, m);

	W = imatrix(n, m);
	R = dmatrix(n, n);
	mid(A, R);
	mid(B, Bm);
	if (inv(R)) {
		Xc = R * Bm; //% midpoint solution
		C = reduce(R * A); //% preconditioning
		Z = reduce(R * (B - A * Xc)); //% preconditioning and residual correction
		rad(Z, Zr); //% zr = rad(z)
		rad(C, Cr); //% Cr = rad(C)
		if (rhoSpectral(Cr) >= 1) {
			return false; //% the necessary condition is not fulfilled
		}
		mag(Xc, Xcm); //% abs of the midpoint solution
		Unit(I);
		M = I - Cr;
		if (inv(M)) { //% inverse of M
			Xast = M * Zr;
			W = Xc + Xast * interval(-1.0, 1.0);
			return true;
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauserSkeel <<<

//%============================================================================
//%
//% Methods for solving interval linear systems
//%
//%============================================================================

bool
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<< 
//% for solving interval linear systems standard version (after Rohn)
//%----------------------------------------------------------------------------
BauerSkeel(const imatrix&	A, //% interval matrix
	const ivector&			b, //% interval right-hand vector
	ivector&				w) //% interval solution vector
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), M(n, n), Cr(n, n), Rm(n, n), Ar(n, n);
	dvector		bm(n), xc(n), zr(n), xcm(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		xc = R * bm; //% tilde x - midpoint solution
		C = R * A; //% preconditioning
		z = R * b; //% preconditioning
		inf(C, M); //% M=inf([A])=I-|R|*rad([A])
		rad(C, Cr); //% Cr=rad([C])=|R|*rad([A])
		if (rhoSpectral(Cr) >= 1) {
			return false; //% the necessary condition is not fulfilled
		}
		rad(z, zr); //% zr=rad([z])=|R|*rad([b])
		mag(xc, xcm); //% xc=|xc|
		if (inv(M)) { //% Rc=M^{-1}
			w = xc + M * (Cr * xcm + zr) * interval(-1.0, 1.0);
			return true; //% success!!!
		}
		return false;
	}
	return false; //% singular matrix
} //% >>> BauerSkeel <<<

bool
//%----------------------------------------------------------------------------
//% >>> BauerSkeel method <<<
//% with preconditioning for solving interval linear systems 
//%----------------------------------------------------------------------------
BauerSkeelP(const imatrix& A, //% interval matrix
	const ivector& b,		 //% interval right-hand vector
	ivector& w)				 //% interval solution vector
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), M(n, n), Cr(n, n);
	dvector		bm(n), xc(n), zr(n), xcm(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) {
		C = R * A; //% preconditioning
		z = R * b; //% preconditioning
		return BauerSkeel(C, z, w);
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% BauerSkeel method for solving interval linear systems
//% standard version (after Rohn). The method takes precalculated 
//% preconditioner
//%----------------------------------------------------------------------------
BauerSkeel(const dmatrix& R, //% real nonsingular matrix - preconditioner
	const imatrix& A,		 //% interval matrix
	const ivector& b,		 //% interval right-hand vector
	ivector& w)				 //% interval solution vector
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	dmatrix		M(n, n), Cr(n, n);
	dvector		bm(n), xc(n), zr(n), xcm(n);
	ivector		z(n);
	imatrix		C(n, n);

	w = ivector(n);
	mid(b, bm); //% midpoint vector
	xc = R * bm; //% tilde x - midpoint solution
	C = R * A; //% preconditioning
	z = R * b; //% preconditioning
	inf(C, M); //% M=inf([A])=I-|R|*rad([A])
	rad(C, Cr); //% Cr=rad(C)=|R|*rad([A])
	if (rhoSpectral(Cr) >= 1) {
		return false; //% the necessary condition is not fulfilled
	}
	rad(z, zr); //% zr=rad([z])=|R|*rad([b])
	mag(xc, xcm); //% xc=|xc|
	if (inv(M)) { //% Rc=M^{-1}
		w = xc + M * (Cr*xcm + zr) * interval(-1.0, 1.0);
		return true; //% success!!!
	}
	return false; //% singular matrix
}

bool
BauerSkeelRC(const imatrix& A, //% interval matrix
	const ivector& b,		   //% interval right-hand vector
	ivector& w)				   //% interval solution
//%----------------------------------------------------------------------------
//% BauerSkeel method for solving interval linear systems
//% with residual correction
//%----------------------------------------------------------------------------
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n);
	dvector		bm(n), xc(n);

	w = ivector(n);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse - preconditioner
		xc = R * bm; //% tilde x - midpoint solution
		z = b - A * xc; //% residual correction
		if (BauerSkeel(R, A, z, w)) {
			w = xc + w;
			return true; //% success!!!
		}
		return false;
	}
	return false; //% singular matrix
}

//%================================================================================
//% Other methods
//%================================================================================

bool
BauerSkeel(const ivector& p,
	const omatrix& oA,
	const ovector& ob,
	ivector& w)
	//%----------------------------------------------------------------------------
	//% Hladik's version of the Bauer-Skeel method for solving parametric 
	//% interval linear systems
	//%----------------------------------------------------------------------------
{
	int	n = ob.size(), np = p.size();
	dvector	pm(np), bm(n), zr(n), x0(n), mx0(n), y(n), w0(n), wup(n), wdown(n);
	dmatrix	R(n, n), M(n, n), Rm(n, n), I(n, n);
	ivector	z(n);
	imatrix	C(n, n);

	mid(p, pm);
	bm = v_from_dep(pm, ob);
	R = m_from_dep(pm, oA);
	if (inv(R)) {
		x0 = R * bm; // x0 = mid(A(p))^{-1}*mib(b)
		c_matrix(R, p, oA, C); // C = mid(A(p))^{-1}*A(p)
		z_vector(x0, R, p, oA, ob, z); // b = mid(A(p))^{-1}*(b(q) - A(p)*x0)
		rad(z, zr); // br = rad(b)
		rad(C, M); // M = rad(C)
		Unit(I);
		Rm = I - M;
		if (inv(Rm)) { // now Rm = (I - M)^{-1}
			y = Rm * zr; // y = (I - M)^{-1}*rad(mid(A(p))^{-1}*(b(q) - A(p)*x0))
			wdown = x0 - y;
			wup = x0 + y;
			for (int i = 0; i < n; i++) {
				w[i] = interval(wdown[i], wup[i]);
			}
			return true;
		}
	}
	return false;
}

bool
HladikImpr(const aafrmatrix& A,
	const aafrvector& b,
	ivector& w)
	//%----------------------------------------------------------------------------
	//% BauerSkeel version with calculation performed using
	//% reduced affine arithmetic
	//% Improved bounds
	//% From Milan's paper "Enclosures..."
	//%----------------------------------------------------------------------------
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), Cm(n, n), M(n, n), I(n, n), Rm(n, n);
	dvector		bm(n), x0(n), zr(n), y(n), wdown(n), wup(n);
	AAFR		s, s2;

	w = ivector(n);
	mid(A, R);
	mid(b, bm);
	if (inv(R)) {
		x0 = R * bm;
		for (int i = 0; i < n; i++) {
			s2 = 0;
			for (int j = 0; j < n; j++) {
				s = 0.0;
				for (int k = 0; k < n; k++)
					s = s + A(j, k) * x0[k];
				s2 = s2 + R(i, j) * (b[j] - s);
			}
			z[i] = s2.reduce();
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				s = 0.0;
				for (int k = 0; k < n; k++) {
					s = s + R(i, k) * A(k, j);
				}
				C(i, j) = s.reduce();
			}
		}
		rad(z, zr); // br = rad(z)
		rad(C, M); // M = rad(C)
		Unit(I);
		Rm = I - M;
		if (inv(Rm)) { // now Rm = (I - M)^{-1}
					   //y = Rm * zr; // y = (I - M)^{-1}*rad(mid(A(p))^{-1}*(b(q) - A(p)*x0))
					   //wdown = x0 - y;
					   //wup = x0 + y;
					   //for (int i = 0; i < n; i++) {
					   //	w[i] = interval(wdown[i], wup[i]);
					   //}
			w = x0 + Rm * zr * interval(-1.0, 1.0);
			return true;
		}
	}
	return false; // singular matrix
}

bool
Hladik(const aafmatrix& A,
	const aafvector& b,
	ivector& w)
	//%----------------------------------------------------------------------------
	//% BauerSkeel version with calculation performed using
	//% affine arithmetic
	//%----------------------------------------------------------------------------
{
	if (!A.is_square()) return false;
	int	n = A.num_rows();
	ivector 	z(n);
	imatrix	C(n, n);
	dmatrix	R(n, n), Cm(n, n), M(n, n), I(n, n), Rm(n, n);
	dvector	bm(n), x0(n), zr(n), y(n), wdown(n), wup(n);
	AAF	s, s2;

	w = ivector(n);
	mid(A, R);
	mid(b, bm);
	if (inv(R)) {
		x0 = R * bm;
		for (int i = 0; i < n; i++) {
			s2 = 0;
			for (int j = 0; j < n; j++) {
				s = 0;
				for (int k = 0; k < n; k++) {
					s = s + A(j, k) * x0[k];
				}
				s2 = s2 + R(i, j) * (b[j] - s);
			}
			z[i] = s2.reduce();
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				s = 0;
				for (int k = 0; k < n; k++)
					s = s + R(i, k) * A(k, j);
				C(i, j) = s.reduce();
			}
		}
		rad(z, zr); // br = rad(z)
		rad(C, M); // M = rad(C)
		Unit(I);
		Rm = I - M;
		if (inv(Rm)) { // now Rm = (I - M)^{-1}
			y = Rm * zr; // y = (I - M)^{-1}*rad(mid(A(p))^{-1}*(b(q) - A(p)*x0))
			wdown = x0 - y;
			wup = x0 + y;
			for (int i = 0; i < n; i++)
				w[i] = interval(wdown[i], wup[i]);
			return true;
		}
	}
	return false; // singular matrix
}

bool
Hladik(const ivector& p,
	const ivector& q,
	const omatrix& oA,
	const ovector& ob,
	ivector& w)
	//%----------------------------------------------------------------------------
	//% Method based on Popova's characterisation of param sol set
	//%----------------------------------------------------------------------------
{
	int			n = ob.size(), np = p.size(), nq = q.size();
	dvector		pm(np), qm(nq), bm(n), zr(n), x0(n), mx0(n), y(n), w0(n), wup(n), wdown(n);
	dmatrix		R(n, n), M(n, n), Rm(n, n), I(n, n);
	ivector		z(n);
	imatrix		C(n, n);

	mid(p, pm);
	mid(q, qm);
	bm = v_from_dep(qm, ob);
	R = m_from_dep(pm, oA);
	if (inv(R)) {
		x0 = R * bm; // x0 = mid(A(p))^{-1}*mib(b)
		c_matrix(R, p, oA, C); // C = mid(A(p))^{-1}*A(p)
		z_vector(x0, R, p, q, oA, ob, z); // b = mid(A(p))^{-1}*(b(q) - A(p)*x0)
		rad(z, zr); // br = rad(b)
		rad(C, M); // M = rad(C)
		Unit(I);
		Rm = I - M;
		if (inv(Rm)) { // now Rm = (I - M)^{-1}
			y = Rm * zr; // y = (I - M)^{-1}*rad(mid(A(p))^{-1}*(b(q) - A(p)*x0))
			wdown = x0 - y;
			wup = x0 + y;
			for (int i = 0; i < n; i++)
				w[i] = interval(wdown[i], wup[i]);
			return true;
		}
	}
	return false;
}

bool
HladikImpr(const ivector& x,
	const ivector& p,
	const omatrix& oA,
	const ovector& ob,
	ivector& w)
{
	int			n = ob.size(), np = p.size();
	dvector		x0(n), yk(n), bk(n), pm(np), bm(n), y(n), z(n), zabs(n), zi(n), wup(n), wdown(n), w0(n);
	dmatrix		R(n, n), Ak(n, n), Y(n, n), Yabs(n, n), Z(n, n), Zabs(n, n), Zi(n, n), M(n, n), I(n, n);
	ivector		ak(n);
	imatrix		C(n, n);

	w = ivector(n);
	mid(p, pm);
	bm = v_from_dep(pm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	Z.fill_in(0);
	Y.fill_in(0);
	y.fill_in(0);
	z.fill_in(0);
	if (inv(R)) {
		x0 = R * bm; // to jest x*
		for (int k = 0; k < np; k++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					Ak(i, j) = oA(i, j)[k];
				}
				bk[i] = ob[i][k];
			}
			ak = R * (Ak * x - bk);
			if (ak >= 0) {
				Y = Y + (R * Ak) * p[k].rad();
				y = y + (R * (Ak * x0 - bk)) * p[k].rad();
			}
			else if (ak <= 0) {
				Y = Y - (R * Ak) * p[k].rad();
				y = y - (R * (Ak * x0 - bk)) * p[k].rad();
			}
			else {
				Zi = R * Ak;
				mag(Zi, Zabs);
				Z = Z + Zabs * p[k].rad();
				zi = R * (Ak * x0 - bk);
				mag(zi, zabs);
				z = z + zabs * p[k].rad();
			}
		}
		mag(Y, Yabs);
		Unit(I);
		M = I - Yabs - Z;
		if (inv(M)) {
			w0 = M * (y + z);
			wdown = x0 - w0;
			wup = x0 + w0;
			for (int i = 0; i < n; i++)
				w[i] = interval(wdown[i], wup[i]);
			return true;
		}
	}
	return false;
}

bool
HladikImpr2(const ivector& x,
	const ivector& p,
	const omatrix& oA,
	const ovector& ob,
	ivector& w)
{
	int			n = ob.size(), np = p.size();
	dvector		x0(n), yk(n), bk(n), pm(np), bm(n), y(n), z(n), zi(n), wup(n), wdown(n), w0(n);
	dmatrix		R(n, n), Ak(n, n), Y(n, n), Yabs(n, n), Z(n, n), Zi(n, n), M(n, n), I(n, n), RAkp(n, n);
	ivector		ak(n);
	imatrix		C(n, n);

	w = ivector(n);
	mid(p, pm);
	bm = v_from_dep(pm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	Z.fill_in(0);
	Y.fill_in(0);
	y.fill_in(0);
	z.fill_in(0);
	if (inv(R)) {
		x0 = R * bm; // to jest x*
		for (int k = 0; k < np; k++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++)
					Ak(i, j) = oA(i, j)[k];
				bk[i] = ob[i][k];
			}
			ak = R * (Ak * x - bk);
			RAkp = (R * Ak) * p[k].rad();
			for (int i = 0; i < n; i++) {
				if (ak[i] >= 0.0) {
					for (int j = 0; j < n; j++)
						Y(i, j) = Y(i, j) + RAkp(i, j);
					y[i] = y[i] + (R * (Ak * x0 - bk))[i] * p[k].rad();
				}
				else if (ak[i] <= 0.0) {
					for (int j = 0; j < n; j++)
						Y(i, j) = Y(i, j) - RAkp(i, j);
					y[i] = y[i] - (R * (Ak * x0 - bk))[i] * p[k].rad();
				}
				else {
					for (int j = 0; j < n; j++) {
						Z(i, j) = Z(i, j) + fabs(RAkp(i, j));
					}
					z[i] = z[i] + fabs((R * (Ak * x0 - bk))[i]) * p[k].rad();
				}
			}
		}
		mag(Y, Yabs);
		Unit(I);
		M = I - Yabs - Z;
		if (inv(M)) {
			w0 = M * (y + z);
			wdown = x0 - w0;
			wup = x0 + w0;
			for (int i = 0; i < n; i++)
				w[i] = interval(wdown[i], wup[i]);
			return true;
		}
	}
	return false;
}