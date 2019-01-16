//%--------------------------------------------------------------
//%
//%  Different versions of Skalna's Direct Method
//%
//%--------------------------------------------------------------
#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../utils/randnum.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../truss/TrussStructure.h"
#include "../param_solver/iterative.h"
#include "mism.h"

bool 
h_matrix(
const dmatrix& A //% revised affine matrix
)
//%--------------------------------------------------------------
//%
//% Check if A is an H-matrix
//% Def: A is an H-matrix if there exists a real vector u > 0,
//% such that <A>u > 0, where <A> - Ostrovsky comparison matrix
//% the matrix A in argument is usually the Ostrovsky matrix
//%
//%--------------------------------------------------------------
{
	if (!A.is_square()) return false;
	int		n = A.num_rows();
	dvector	e(n), x(n);

	e.fill_in(1);
	if (gw_gausse(A, e, x)) { //% rozwiazanie ukladu <A>x = e
		for (int i = 0; i < n; i++) {
			if (x[i] <= 0.0) {
				return false; //% not decided
			}
		}
		return true; //% is H-matrix
	}
	return false; //% not H-matrix (singularity)
}

bool 
ISM(
const aafrmatrix&	A, //% revised affine matrix
const aafrvector&	b, //% revised right-hand vector
ivector&			w  //% interval solution
)
//%--------------------------------------------------------
//%
//% Direct Method
//% x^OI = xc + <RA>^{-1}|R(b-A*xc)|[-1,1]
//%
//%--------------------------------------------------------
{
	if (!A.is_square()) return false;
	int		n = A.num_rows();
	ivector	z(n);
	imatrix	B(n, n);
	dmatrix	R(n, n), Bm(n, n);
	dvector	bm(n), xc(n), magz(n), cmz(n);
	
	mid(A, R);  //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		xc = R * bm; //% midpoint solution
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				AAFR si = 0.0;
				for (int k = 0; k < n; k++) {
					if ((A(j, k).reduce() != 0.0) && (xc[k] != 0.0)) {
						si = si + A(j, k) * xc[k];
					}
				}
				s = s + R(i, j) * (b[j] - si);
			}
			z[i] = s.reduce();
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if ((A(k, j).reduce() != 0.0) && (R(i, k) != 0.0)) {
						s = s + R(i, k) * A(k, j);
					}
				}
				B(i, j) = s.reduce(); //% R*A
			}
		}
		mig(B, Bm);
		if (!h_matrix(Bm)) { //% H-matrix property verification
			return false;
		}
		mag(z, magz); //% Ostrowski comparison matrix
		if (gw_gausse(Bm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				w[i] = xc[i] + cmz[i] * interval(-1.0, 1.0);
			}
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

bool
ISM(
const aafrmatrix&	A, //% revised affine matrix
const aafrmatrix&	B, //% revised right-hand matrix
imatrix&			w  //% interval solution
)
//%--------------------------------------------------------
//%
//% Direct Method
//% for multiple right-hand side (MRHS)
//% X^OI = XC + <RA>^{-1}|R(B-A*XC)|[-1,1]
//%
//%--------------------------------------------------------
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	int			m = B.num_cols();
	imatrix		z(n, m), C(n, n);
	aafrmatrix	Ax(n, m);
	dmatrix		R(n, n), Cm(n, n);
	dmatrix		bm(n,m), x(n,m), magz(n,m), cmz(n,m);

	w = imatrix(n,m);
	w.fill_in(0);
	mid(A, R); //% midpoint matrix
	mid(B, bm); //% midpoint right-hand matrix
	if (inv(R)) { //% midpoint inverse (must be non-singular)
		x = R * bm; //% midpoint solution
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if ((A(i, k).reduce() != 0.0) && (x(k, j) != 0.0)) {
						s = s + A(i, k) * x(k, j);
					}
				}
				Ax(i, j) = B(i, j) - s; //% residual correction
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (R(i, k) != 0.0 && Ax(k, j).reduce() != 0.0) {
						s = s + R(i, k) * Ax(k, j);
					}
				}
				z(i, j) = s.reduce();
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (R(i, k) != 0.0 && Ax(k, j).reduce() != 0.0) {
						s = s + R(i, k) * A(k, j);
					}
				}
				C(i, j) = s.reduce();
			}
		}
		mig(C, Cm);
		if (!h_matrix(Cm)) { //% H-matrix property verification
			return false; //% not an H-matrix
		}
		mag(z, magz);
		if (inv(Cm)) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					real ss = 0.0;
					for (int k = 0; k < n; k++) {
						ss = ss + Cm(i, k) * magz(k, j);
					}
					w(i, j) = x(i, j) + ss * interval(-1.0, 1.0);
				}
			}
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

bool
ISM(
const aafrmatrix&	A, //% revised affine matrix
const aafrmatrix&	B, //% revised affine right-hand matrix
imatrix&			w, //% interval solution
dmatrix&			R  //% midpoint inverse
)
//%--------------------------------------------------------
//%
//% Direct Method
//% for multiple right-hand side (MRHS)
//% X^OI = XC + <RA>^{-1}|R(B-A*XC)|[-1,1]
//% Additionally midpoint inverse is returned
//%
//%--------------------------------------------------------
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	int			m = B.num_cols();
	imatrix		z(n, m), C(n, n);
	aafrmatrix	Ax(n, m);
	dmatrix		Cm(n, n);
	dmatrix		bm(n, m), x(n, m), magz(n, m), cmz(n, m);

	w = imatrix(n, m);
	w.fill_in(0);
	R = dmatrix(n, n);
	mid(A, R); //% midpoint matrix
	mid(B, bm); //% midpoint right-hand matrix
	if (inv(R)) { //% midpoint inverse
		x = R * bm; // midpoint solution
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (A(i, k).reduce() != 0.0 && x(k, j) != 0.0) {
						s = s + A(i, k) * x(k, j);
					}
				}
				Ax(i, j) = B(i, j) - s;
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (R(i, k) != 0.0 && Ax(k, j).reduce() != 0.0) {
						s = s + R(i, k) * Ax(k, j);
					}
				}
				z(i, j) = s.reduce();
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (R(i, k) != 0.0 && A(k, j).reduce() != 0.0) {
						s = s + R(i, k) * A(k, j);
					}
				}
				C(i, j) = s.reduce();
			}
		}
		mig(C, Cm);
		if (!h_matrix(Cm)) { //% H-matrix property verification
			return false; //% not an H-matrix
		}
		mag(z, magz);
		if (inv(Cm)) {
			//if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					real ss = 0.0;
					for (int k = 0; k < n; k++) {
						ss = ss + Cm(i, k) * magz(k, j);
					}
					w(i, j) = x(i, j) + ss * interval(-1.0, 1.0);
				}
			}
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

bool 
ISM(
const aafrmatrix&	A, //% revised affine matrix
const aafrvector&	b, //% revised right-hand vector
ivector&			w, //% p-solution (revised affine form)
dmatrix&			R  //% midpoint inverse
)
//%--------------------------------------------------------
//%
//% Direct Method
//% Additionally midpoint inverse is returned
//%
//%--------------------------------------------------------
{
	int		n = A.num_rows();
	ivector	z(n);
	imatrix	C(n, n);
	dmatrix	Cm(n, n);
	dvector	bm(n), x(n), magz(n), cmz(n);

	assert(A.is_square());

	w.fill_in(0.0);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if(inv(R)) { //% midpoint inverse
		x = R * bm; //% midpoint solution
		for(int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for(int j = 0; j < n; j++) {
				AAFR si = 0.0;
				for(int k = 0; k < n; k++) {
					if ((A(j, k).reduce() != 0.0) && (x[k] != 0.0)) {
						si = si + A(j, k) * x[k];
					}
				}
				s = s + R(i, j) * (b[j] - si);
			}
			z[i] = s.reduce();
		}
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for(int k = 0; k < n; k++) {
					if ((R(i, k) != 0.0) && (A(k, j).reduce() != 0.0)) {
						s = s + R(i, k) * A(k, j);
					}
				}
				C(i, j) = s.reduce();
			}
		}
		mig(C, Cm);
		if(!h_matrix(Cm)) { //% H-matrix property verification
			return false; //% not an H-matrix
		}
		mag(z, magz);
		if(gw_gausse(Cm, magz, cmz)) {
			for(int i = 0; i < n; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0);
			}
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

bool
ISMP0(
const aafrmatrix&	A, //% revised affine matrix
const aafrvector&	b, //% revised right-hand vector
ivector&			w, //% p-solution (revised affine form)
dmatrix&			R  //% midpoint inverse
)
//%--------------------------------------------------------
//%
//% Direct Method with post-conditioning
//% Additionally midpoint inverse is returned
//%
//%--------------------------------------------------------
{
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		Cm(n, n);
	dvector		bm(n), x(n), magz(n), cmz(n);
	aafrmatrix	B(n, n);

	assert(A.is_square());

	w.fill_in(0.0);
	mid(A, R); //% revised affine matrix
	mid(b, bm); //% revised right-hand vector
	if (inv(R)) { //% midpoint inverse
		B = A * R; //% post-conditioning
		x = bm; //% midpoint solution
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				AAFR si = 0.0;
				for (int k = 0; k < n; k++) {
					if ((B(j, k) == 0.0) || (x[k] == 0.0)) {
						continue;
					}
					si = si + B(j, k) * x[k];
				}
				s = s + b[j] - si; //% residual correction
			}
			z[i] = s.reduce(); //% right-hand vector with residual correction
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C(i, j) = B(i, j).reduce();
			}
		}
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		if (!h_matrix(Cm)) { //% if not an H-matrix, the method does not apply
			return false;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0);
			}
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

bool
ISMP(
const aafrmatrix&	A, //% revised affine matrix
const aafrvector&	b, //% revised right-hand vector
ivector&			w, //% p-solution (revised affine form)
dmatrix&			R  //% midpoint inverse
)
//%--------------------------------------------------------
//%
//% Direct Method with post-conditioning
//% Additionally midpoint inverse is returned
//%
//%--------------------------------------------------------
{
	int			n = A.num_rows();
	ivector		z(n);
	imatrix		C(n, n);
	dmatrix		Cm(n, n);
	dvector		bm(n), x(n), magz(n), cmz(n);
	aafrmatrix	B(n, n);

	assert(A.is_square());

	w.fill_in(0.0);
	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		B = A * R; //% post-conditioning
		x = bm; //% midpoint solution
		for (int i = 0; i < n; i++) {
			AAFR s = 0.0;
			for (int j = 0; j < n; j++) {
				AAFR si = 0.0;
				for (int k = 0; k < n; k++) {
					if ((B(j, k).reduce() != 0.0) && (x[k] != 0.0)) {
						si = si + B(j, k) * x[k];
					}
				}
				s = s + b[j] - si; //% residual correction
			}
			z[i] = s.reduce(); //% right-hand vector with residual correction
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				C(i, j) = B(i, j).reduce();
			}
		}
		mig(C, Cm); //% Ostrovski comparison matrix for C=B*R
		if (!h_matrix(Cm)) { //% if not an H-matrix, the method does not apply
			return false;
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) {
			for (int i = 0; i < n; i++) {
				w[i] = x[i] + cmz[i] * interval(-1.0, 1.0);
			}
			w = R * w;
			return true;
		}
	}
	return false; //% singular midpoint matrix => A singular
}

void 
OUSystem(
const aafrmatrix&	A,  //% revised affine matrix
const aafrvector&	b,  //% revised affine right-hand vector
aafrmatrix&			AA, //% augmented revised affine matrix 
aafrvector&			bb  //% augmented revised affine right-hand vector
) 
//%--------------------------------------------------------------
//%
//% Create augmented system 
//% from over- or under-determined system
//%
//%--------------------------------------------------------------
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
}
