//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Generalized Expansion Method (GEM)
//%     for solving parametric interval linear systems
//%     GEM is based on:
//%     a) computing inverse of an affine matrix by using Neumann series
//%     b) computing inverse of an affine matrix using Rohn's formula
//%
//%##########################################################################
//%

#include "../utils/stdafx.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../param_solver/direct_solvers.h"
#include "../truss/TrussStructure.h"
#include "../param_solver/kolev.h"
#include "../examples/pexamples.h"

int
//%----------------------------------------------------------------------------
//% Generalized Expansion Method (GEM), basic version based on computing 
//% inverse of an affine matrix by using Neumann series
//%----------------------------------------------------------------------------
gem(int m, //% order of the method
	int op, //% option: 0 - with residual correction, 1 - without residual correction
	const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& w) //% revised affine solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix C(n, n), Cinv(n, n), mA(n, n), SA(n, n), mAn(n, n);
	aafrvector z(n);
	dmatrix R(n, n), Cr(n, n), I(n, n), D(n, n), Zn(n, n);
	dvector xc(n), bc(n);
	imatrix Ci(n, n);

	w = aafrvector(n);
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		mid(b, bc); //% midpoint vector
		xc = R * bc; //% midpoint solution
		if (op == 0) {
			z = R * (b - A * xc); //% left pre-conditioning and residual correction
		}
		else {
			z = R * b; //% left pre-conditioning
		}
		C = R * A; //% left preconditioning
		Ci = reduce(C);
		rad(Ci, Cr); //% radius of C
		Unit(I);
		D = I - Cr; //% (I-C^delta)
		if (inv(D)) { //% (I-C^delta)^{-1}
			Zn = Cr;
			//% this can be improved to have lover complexity
			for (int j = 0; j < m; j++) {
				Zn = Zn * Cr;
			}
			Zn = Zn * D; //% Zn is OK
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					mA(i, j) = I(i, j) - C(i, j); //% I-C
				}
			}
			mAn = mA;
			SA.fill_in(0.0);
			//% here we can gain some time
			for (int j = 0; j < m; j++) {
				SA = SA + mAn;
				mAn = mAn * mA;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					real zn = Zn(i, j);
					Cinv(i, j) = I(i, j) + SA(i, j);
					Cinv(i, j) = Cinv(i, j) + interval(-zn, zn);
				}
			}
			w = Cinv * z; //% final solution without residual correction
			if (op == 0) {
				for (int i = 0; i < n; i++) {
					w[i] = xc[i] + w[i]; //% for residual correction
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Generalized Expansion Method, basic version based on computing inverse
//% of an affine matrix by using Neumann series. Additionally the midpoint
//% inverse is returned in argument R.
//%----------------------------------------------------------------------------
gem(int m, //% order of the method
	int op, //% option: 0 - with residual correction, 1 - without residual correction
	const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand vector
	aafrvector& w, //% revised affine solution vector
	dmatrix& R) //% midpoint inverse
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix C(n, n), Cinv(n, n), AA(n, n), mA(n, n), SA(n, n), mAn(n, n);
	aafrvector bb(n), z(n), zf(n);
	dmatrix Z(n, n), Cr(n, n), I(n, n), D(n, n), Zn(n, n);
	dvector xc(n), bc(n);
	imatrix Ci(n, n);

	mid(b, bc); //% midpoint vector
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		xc = R * bc;
		if (op == 0) {
			z = R * (b - A * xc); //% left pre-conditioning and residual correction
		}
		else {
			z = R * b; //% left pre-conditioning
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					s = s + R(i, k) * A(k, j);
				}
				C(i, j) = s; //% R * A
			}
		}
		Unit(I);
		Ci = reduce(C);
		rad(Ci, Cr); //% radius of C
		D = I - Cr; //% (I-C^delta)
		if (inv(D)) { //% (I-C^delta)^{-1}
			Zn = Cr;
			for (int j = 0; j < m; j++) {
				Zn = Zn * Cr;
			}
			Zn = Zn * D; //% Zn is OK
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					mA(i, j) = I(i, j) - C(i, j);
				}
			}
			mAn = mA;
			SA.fill_in(0.0);
			for (int j = 0; j < m; j++) {
				SA = SA + mAn;
				mAn = mAn * mA;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					real zn = Zn(i, j);
					Cinv(i, j) = I(i, j) + SA(i, j);
					Cinv(i, j) = Cinv(i, j) + interval(-zn, zn);
				}
			}
			w = Cinv * z; //% final solution
			if (op == 0) {
				for (int i = 0; i < n; i++) {
					w[i] = xc[i] + w[i];
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Generalized Expansion method, basic version based on computing inverse
//% of an affine matrix by using Neumann series. Method solves systems with
//% multiple right-hand side.
//%----------------------------------------------------------------------------
gem(int mm, //% order of the method
	int op, //% option: 0 - with residual correction, 1 - without residual correction
	const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& b, //% revised right-hand matrix
	aafrmatrix& w) //% revised affine solution matrix
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = b.num_cols();
	aafrmatrix C(n, n), Cinv(n, n), mA(n, n), SA(n, n), mAn(n, n);
	aafrmatrix z(n, m), Ax(n, m);
	imatrix CC(n, n);
	dmatrix R(n, n), Z(n, n), Cr(n, n), I(n, n), D(n, n), Zn(n, n);
	dmatrix xc(n, m), bc(n, m);

	w = aafrmatrix(n, m);
	mid(b, bc); //% midpoint vector
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		xc = R * bc; //% midpoint solution
		if (op == 0) { //% left pre-conditioning and resudial correction
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					AAFR s = 0.0;
					for (int k = 0; k < n; k++) {
						if (A(i, k) == 0.0) continue;
						s = s + A(i, k) * xc(k, j);
					}
					Ax(i, j) = b(i, j) - s;
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					AAFR s = 0.0;
					for (int k = 0; k < n; k++) {
						s = s + R(i, k) * Ax(k, j);
					}
					z(i, j) = s;
				}
			}
		}
		else {
			z = R * b; //% left pre-conditioning
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					s = s + R(i, k) * A(k, j);
				}
				C(i, j) = s; //% R * A
				CC(i, j) = s.reduce(); //% [R*A]
			}
		}
		Unit(I);
		rad(CC, Cr); //% radius of C
		D = I - Cr; //% (I-C^delta)
		if (inv(D)) { //% (I-C^delta)^{-1}
			Zn = Cr;
			for (int j = 0; j < mm; j++) {
				Zn = Zn * Cr;
			}
			Zn = Zn * D;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					mA(i, j) = I(i, j) - C(i, j);
				}
			}
			mAn = mA;
			SA.fill_in(0.0);
			for (int j = 0; j < mm; j++) {
				SA = SA + mAn;
				mAn = mAn * mA;
			}

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					real zn = Zn(i, j);
					Cinv(i, j) = I(i, j) + SA(i, j);
					Cinv(i, j) = Cinv(i, j) + interval(-zn, zn);
				}
			}
			w = Cinv * z; //% final solution
			if (op == 0) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						w(i, j) = xc(i, j) + w(i, j);
					}
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Generalized Expansion method based on computing inverse of an affine 
//% matrix by using Rohn's formula. The method has lower time complexity
//% than the basic version.
//%----------------------------------------------------------------------------
gemMH(int m, //% order of the method
	const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix C(n, n), Cinv(n, n), V0(n, n), mA(n, n), SA(n, n), mAn(n, n);
	aafrvector z(n), zf(n), mZ(n), smZ(n), smZ2(n);
	imatrix CC(n, n), H(n, n);
	dmatrix R(n, n), Z(n, n), Cr(n, n), I(n, n), M(n, n), Zn(n, n);
	dvector xc(n), bc(n);

	mid(b, bc); //% midpoint vector
	mid(A, R); //% midpoint matrix
	if (inv(R)) { //% midpoint inverse
		xc = R * bc;
		z = R * (b - A * xc); //% left pre-conditioning and residual correction
		C = R * A; //% left pre-conditioning
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V0(i, j) = I(i, j) - C(i, j);
				CC(i, j) = C(i, j).reduce();
			}
		}
		rad(CC, Cr); //% radius of C
		M = I - Cr; //% (I-C^delta) 
		if (inv(M)) { //% (I-C^delta)^{-1} >= 0
			for (int i = 0; i < n; i++) {  //% Computing H matrix, H=inv(C)
				for (int j = 0; j < n; j++) {
					if (i == j) {
						H(i, j) = interval(M(i, j) / (2.0 * M(i, j) - 1.0), M(i, j));
					}
					else {
						H(i, j) = interval(-M(i, j), M(i, j));
					}
				}
			}
			mZ = V0 * z;
			smZ.fill_in(0.0); //% V0^m * z
			for (int j = 0; j < m; j++) {
				smZ = smZ + mZ;
				mZ = V0 * mZ;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {					
					real hr = H(i, j).rad();
					Cinv(i, j) = H(i, j).mid();
					Cinv(i, j) = Cinv(i, j) + interval(-hr, hr);
				}
			}
			w = z + smZ + Cinv * mZ; //% final solution
			for (int i = 0; i < n; i++) {
				w[i] = xc[i] + w[i];
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Generalized Expansion method based on computing inverse of an affine 
//% matrix by using Rohn's formula. The method has lower time complexity than
//% the basic version. The method solves systems with multiple rihgt-hand side.
//%----------------------------------------------------------------------------
gemMH(int mm, //% order of the method
	const aafrmatrix& A, //% revised affine matrix
	const aafrmatrix& b, //% revised affine right-hand side matrix
	aafrmatrix& w) //% interval solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows(), m = b.num_cols();
	aafrmatrix C(n, n), Cinv(n, n), V0(n, n), mA(n, n);
	aafrmatrix z(n, m), mZ(n, m), smZ(n, m), smZ2(n, m), Ax(n, m);
	imatrix CC(n, n), H(n, n);
	dmatrix R(n, n), Z(n, n), Cr(n, n), I(n, n), M(n, n), Zn(n, n);
	dmatrix xc(n, m), bc(n, m);

	mid(A, R);
	if (inv(R)) {
		mid(b, bc);
		xc = R * bc;
		//% residual correction
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					if (A(i, k) == 0.0) continue;
					s = s + A(i, k) * xc(k, j);
				}
				Ax(i, j) = b(i, j) - s;
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				AAFR s = 0.0;
				for (int k = 0; k < n; k++) {
					s = s + R(i, k) * Ax(k, j);
				}
				z(i, j) = s;
			}
		}
		C = R * A; //% left pre-conditioning
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V0(i, j) = I(i, j) - C(i, j);
			}
		}
		CC = reduce(C);
		rad(CC, Cr); //% radius of C
		M = I - Cr; //% (I-C^delta) 
		if (inv(M)) { //% (I-C^delta)^{-1} >= 0
			for (int i = 0; i < n; i++) { //% Computing H matrix, H=inv(C) using Rohn's formula
				for (int j = 0; j < n; j++) {
					if (i == j) {
						H(i, j) = interval(M(i, j) / (2.0 * M(i, j) - 1.0), M(i, j));
					}
					else {
						H(i, j) = interval(-M(i, j), M(i, j));
					}
				}
			}
			mZ = V0 * z;
			smZ.fill_in(0.0); //% V0^m * z
			for (int j = 0; j < mm; j++) {
				smZ = smZ + mZ;
				mZ = V0 * mZ;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					real hr = H(i, j).rad();
					Cinv(i, j) = H(i, j).mid();
					Cinv(i, j) = Cinv(i, j) + interval(-hr, hr);
				}
			}
			w = z + smZ + Cinv * mZ; //% final solution
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					w(i, j) = xc(i, j) + w(i, j);
				}
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

int
//%----------------------------------------------------------------------------
//% Generalized Expansion method, based on Neumann series. The method is
//% a modification of the basic version. The modification aims to lower time
//% complexity.
//%----------------------------------------------------------------------------
gem(int m, //% order of the method
	const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine vector
	aafrvector& w) //% revised affine solution vector
{
	if (!A.is_square()) return NonSquareMatrix;

	int n = A.num_rows();
	aafrmatrix C(n, n), Cinv(n, n), mA(n, n), SA(n, n), mAn(n, n);
	aafrvector z(n), mZ(n), smZ(n);
	imatrix CC(n, n);
	dmatrix R(n, n), Z(n, n), Cr(n, n), I(n, n), D(n, n), Zn(n, n);
	dvector xc(n), bc(n);

	w = aafrvector(n);
	mid(b, bc);
	mid(A, R);
	if (inv(R)) {
		xc = R * bc;
		z = R * (b - A * xc); //% left pre-conditioning and residual correction
		C = R * A; //% left pre-conditioning
		CC = reduce(C);
		Unit(I);
		rad(CC, Cr); //% radius of C
		D = I - Cr; //% (I-C^delta)
		if (inv(D)) { //% (I-C^delta)^{-1}
			Zn = Cr;
			for (int j = 0; j < m; j++) {
				Zn = Zn * Cr;
			}
			Zn = Zn * D;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					mA(i, j) = I(i, j) - C(i, j);
				}
			}
			mZ = mA * z;
			smZ.fill_in(0.0);
			for (int j = 0; j < m; j++) {
				smZ = smZ + mZ;
				mZ = mA * mZ;
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					real zn = Zn(i, j);
					Cinv(i, j) = I(i, j);
					Cinv(i, j) = Cinv(i,j) + interval(-zn, zn);
				}
			}
			w = smZ + Cinv * z; //% final solution
			for (int i = 0; i < n; i++) {
				w[i] = xc[i] + w[i];
			}
			return Success;
		}
		return NotStronglyRegular;
	}
	return SingularMatrix;
}

//%============================================================================
//% Below is the original version of Degrauwe method. There are also 
//% some auxiliary functions
//%============================================================================

void
//%----------------------------------------------------------------------------
//% Converts the matrix oA which is used to compute A,
//% where A_{ij}=oA(i,j)^T*p
//% into the set of matrices which are used to compute 
//% A=\sum_{k=1}^KA^{(k)}p_k
//%----------------------------------------------------------------------------
ConvPtoEps(const omatrix& oA, //% matrix of vectors
	const ivector& p, //% interval parameter
	omatrix& oAn) //% new matrix of vectors
{
	int	n = oA.num_rows(), np = p.size();
	real s;
	dvector	el(np);

	oAn = omatrix(n, n);
	oAn.fill_in(el);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < np; k++) {
				if (k == 0) {
					s = 0.0;
					for (int kk = 1; kk < np; kk++) {
						s += oA(i, j)[kk] * p[kk].mid();
					}
					oAn(i, j)[k] = oA(i, j)[k] + s;
					continue;
				}
				oAn(i, j)[k] = oA(i, j)[k] * p[k].rad();
			}
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Converts the vector ob which is used to compute b,
//% where b_i=ob[i]^T*p
//% into the set of matrices which are used to compute 
//% b=\sum_{k=1}^Kb^{(k)}p_k
//%----------------------------------------------------------------------------
ConvPtoEps(const ovector& ob, //% vector of vectors
	const ivector& q, //% interval vector
	ovector& obn) //% new vector of vectors
{
	int	n = ob.size(), nq = q.size();
	real s;
	dvector el(nq);

	obn = ovector(n);
	obn.fill_in(el);
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < nq; k++) {
			if (k == 0) {
				s = 0.0;
				for (int kk = 1; kk < nq; kk++) {
					s += ob[i][kk] * q[kk].mid();
				}
				obn[i][k] = ob[i][k] + s;
				continue;
			}
			obn[i][k] = ob[i][k] * q[k].rad();
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Computes R*A(p) taking into account dependencies
//%----------------------------------------------------------------------------
RtimesDep(omatrix& oA, //% matrix of vectors
	const dmatrix& R) //% pre-conditioning matrix
{
	int	n = R.num_rows();
	int	np = oA(0, 0).size();
	dmatrix A(n, n);

	for (int k = 0; k < np; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A(i, j) = oA(i, j)[k];
			}
		}
		A = R * A;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				oA(i, j)[k] = A(i, j);
			}
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Computes R*b(p) taking into account dependencies
//%----------------------------------------------------------------------------
RtimesDep(ovector& ob,
	const dmatrix& R)
{
	int	n = R.num_rows();
	int	nq = ob[0].size();
	dvector b(n);

	for (int k = 0; k < nq; k++) {
		for (int i = 0; i < n; i++) {
			b[i] = ob[i][k];
		}
		b = R * b;
		for (int i = 0; i < n; i++) {
			ob[i][k] = b[i];
		}
	}
}

void
//%----------------------------------------------------------------------------
//% COmputes midpoint of a matrix of vectors
//%----------------------------------------------------------------------------
A_mid(const omatrix& oA, const ivector& p, dmatrix& Am)
{
	int	n = oA.num_rows();
	int	np = p.size();
	real s;

	Am = dmatrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			s = 0;
			for (int k = 1; k < np; k++) {
				s += oA(i, j)[k] * p[k].mid();
			}
			Am(i, j) = oA(i, j)[0] + s;
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Converts an element of a matrix of vectors to a real vector
//%----------------------------------------------------------------------------
A_i_oA(int m, const omatrix& oA, dmatrix& Ai)
{
	int n = oA.num_rows();

	Ai = dmatrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Ai(i, j) = oA(i, j)[m];
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Converts an element of a vector of vectors to a real vector
//%----------------------------------------------------------------------------
b_i_ob(int r, const ovector& ob, dvector& bi)
{
	int nq = ob[0].size();

	bi = dvector(nq);
	for (int i = 0; i < nq; i++) {
		bi[i] = ob[r][i];
	}
}

AAFR
//%----------------------------------------------------------------------------
//% Converts real vector to a revised affine form
//%----------------------------------------------------------------------------
vectorToAAFR(const dvector& x)
{
	AAFR a;
	int	n = x.size();
	int	c = 0;

	for (int i = 1; i < n - 1; i++) {
		if (fabs(x[i]) > 0) c++;
	}
	cvector idx(c + 1);
	dvector cfs(c + 1);
	idx[0] = 0; cfs[0] = x[0];
	for (int i = 1, k = 1; i < n - 1; i++) {
		if (fabs(x[i]) > 0) {
			idx[k] = i;
			cfs[k] = x[i];
			k++;
		}
	}
	a = AAFR(idx, cfs, x[n - 1]);
	return a;
}

AAFR
//%----------------------------------------------------------------------------
//% Converts real vector to a revised affine form
//%----------------------------------------------------------------------------
vectorToAAFR(const int np,
	const dvector& x)
{
	AAFR a;
	int	n = x.size();
	int	c = 0;

	for (int i = 0; i < n; i++) {
		if (fabs(x[i]) > 0) c++;
	}
	if (c == 0) {
		return 0.0;
	}
	cvector idx(c);
	dvector cfs(c);
	idx[0] = 0; cfs[0] = x[0];
	int k = 1;
	for (int i = 1; i < n; i++) {
		if (fabs(x[i]) > 0) {
			idx[k] = i + np - 1;
			cfs[k] = x[i];
			k++;
		}
	}
	a = AAFR(idx, cfs);
	return a;
}

void
//%----------------------------------------------------------------------------
//% Verifies 1st oreder Degrauwe's method for solving problems
//%----------------------------------------------------------------------------
main_Degrauwe(Truss* ptruss,
	ivector& w,
	real& t)
{
	clock_t	start, stop;
	AAFR af;
	ivector	p, b, pnew, qnew;
	dvector	el, el2, bi, bii, x0, b0;
	omatrix	oA, oA2, oAn, oB;
	ovector	ob, ob2, obn;
	dmatrix	A, R, Ai, Bi, Ci, C, C2, C3, D2;
	aafrmatrix Bf;
	aafrvector bf;
	imatrix	D, Atmp;

	GlobalStiffnessMatrix(ptruss, A, p, b); //% creates respective vectors and matrices
	oA = OmegaMatrix(A, p, pnew); //% matrix of left hand dependencies
	ob = OmegaVector(b, qnew); //% vector of right hand dependencies

	Atmp = MatrixFromDep(pnew, oA);
	oA2 = oA;
	ob2 = ob;

	ConvPtoEps(oA2, pnew, oAn);
	ConvPtoEps(ob2, qnew, obn);

	ivector pp(pnew.size());

	int n = oA.num_rows();
	A_mid(oAn, pp, R);
	inv(R);

	RtimesDep(oAn, R);
	RtimesDep(obn, R);

	int np = pnew.size(), nq = qnew.size(); //n=oA.num_rows(), 

	oB = omatrix(n, n);
	el = dvector(np + 1);

	C = dmatrix(n, n);
	Ci = dmatrix(n, n);

	Bf = aafrmatrix(n, n);
	bf = aafrvector(n);

	start = clock();

	pp.fill_in(interval(-1, 1));
	pp[0] = 1.0;
	el.fill_in(0.0);
	oB.fill_in(el);
	C.fill_in(0.0);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			oB(i, j)[0] = oAn(i, j)[0]; //B0=I
		}
	}
	//% The estimation of the inverse has the form: B0+sum_{i=1}^KB_i*e_i+B_e*e_e
	//% Now it is assumed that B_ee_e is zero
	//% B0=I, B_i=-A_i, B_e\subseteq (Z^2*(I-Z)^{-1})e_e
	for (int k = 1; k < np; k++) {
		A_i_oA(k, oAn, Ai);
		Bi = Ai * (-1.0);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				oB(i, j)[k] = Bi(i, j);
			}
		}
		mag(Bi, Ci);
		C = C + Ci; // C=sum_{i=1}^np |R*Ai|, C=Z (Z is used in Degrauwes paper)
	}
	//% In this approach A_e=0, A_e normally corresponds to the last element in reduced affine forms
	C2 = C * C; // C^2
	C3 = dmatrix(n, n);
	Unit(C3); // C3=I
	C3 = C3 - C; // C3=I-C
	inv(C3); // (I-C)^{-1}
	D2 = C2 * C3; // B_e=C^2*(I-C)^{-1}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			oB(i, j)[np] = D2(i, j); //Be
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			af = vectorToAAFR(oB(i, j));
			Bf(i, j) = af;
		}
	}
	for (int i = 0; i < n; i++) {
		b_i_ob(i, obn, bi);
		af = vectorToAAFR(np, bi);
		bf[i] = af;
	}
	aafrvector wf = Bf * bf;
	for (int i = 0; i < n; i++) {
		w[i] = wf[i].reduce();
	}
	stop = clock();
	t = difftime(stop, start) / CLOCKS_PER_SEC;
}

void
//%----------------------------------------------------------------------------
//% Verifies 2nd order Degrauwe's method for solving problems
//%----------------------------------------------------------------------------
main_Degrauwe2Order(Truss* ptruss,
	ivector& w,
	real& t)
{
	time_t start;
	AAFR af;
	ivector	p, b, pnew, qnew;
	dvector	el, el2, bi, bii, x0, b0;
	omatrix	oA, oA2, oAn, oB;
	ovector	ob, ob2, obn;
	dmatrix	A, A2, A3, A4, R, Ai, Aii, Bi, Ci, C, C2, C3, D2, A22, Fi, Fii;
	aafrmatrix Bf;
	aafrvector bf;
	imatrix D, Atmp;

	GlobalStiffnessMatrix(ptruss, A, p, b); // creates respective vectors and matrices
	oA = OmegaMatrix(A, p, pnew); // matrix of left hand dependencies
	ob = OmegaVector(b, qnew); // vector of right hand dependencies

	Atmp = MatrixFromDep(pnew, oA);
	oA2 = oA;
	ob2 = ob;

	ConvPtoEps(oA2, pnew, oAn);
	ConvPtoEps(ob2, qnew, obn);

	ivector pp(pnew.size());
	pp.fill_in(interval(-1, 1));
	pp[0] = 1.0;

	int n = oA.num_rows();
	A_mid(oAn, pp, R);
	inv(R);

	RtimesDep(oAn, R);
	RtimesDep(obn, R);

	start = clock();
	int np = pnew.size(), nq = qnew.size(); //n=oA.num_rows(), 

	oB = omatrix(n, n);
	el = dvector(np + 1);
	el.fill_in(0.0);
	oB.fill_in(el);
	C = dmatrix(n, n);
	Ci = dmatrix(n, n);
	C.fill_in(0.0);
	Bf = aafrmatrix(n, n);
	bf = aafrvector(n);
	A2 = dmatrix(n, n);
	A2.fill_in(0.0);
	A3 = dmatrix(n, n);
	A3.fill_in(0.0);
	A4 = dmatrix(n, n);
	A4.fill_in(0.0);
	Fi = dmatrix(n, n);
	Fii = dmatrix(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			oB(i, j)[0] = oAn(i, j)[0]; //B0=I
		}
	}
	for (int k = 1; k < np; k++) {
		A_i_oA(k, oAn, Ai); // extracting Ai
		Bi = Ai * (-1.0); // computing Bi=-Ai
		A22 = Ai * Ai; // A22=Ai^2
		A2 = A2 + A22; // A2=sum of Ai*Ai (i=1,...,n)
		mag(A22, Fi); // Fi=|Ai|
		A3 = A3 + Fi; // A3=sum of |Ai*Ai| (i=1,...,n)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				oB(i, j)[k] = Bi(i, j);
			}
		}
		mag(Bi, Ci);
		C = C + Ci; // C=sum_{i=1}^np |R*Ai| (C=Z from article)
	}
	A2 = A2 * 0.5;
	A3 = A3 * 0.5;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			oB(i, j)[0] = oB(i, j)[0] + A2(i, j); //B0=I+(1/2)sum(Ai^2)
		}
	}
	for (int k = 1; k < np; k++) { // computing 
		for (int kk = k + 1; kk < np; kk++) {
			A_i_oA(k, oAn, Ai);
			A_i_oA(kk, oAn, Aii);
			mag(Ai*Aii + Aii * Ai, Fii); // Fii=|AiAj+AjAi|
			A4 = A4 + Fii; // A4=sum(|AiAj+AjAi|) (i=1,...,n; j=i+1,...,n)
		}
	}
	C2 = C * C*C; // C2=C^3
	C3 = dmatrix(n, n);
	Unit(C3); // C3=I
	C3 = C3 - C; // C3=I-C
	inv(C3); // (I-C)^{-1}
	D2 = A3 + A4 + C2 * C3; // (1/2)sum(|Ai^2|+sum(|AiAj+AjAi|)+C^3*(I-C)^{-1}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			oB(i, j)[np] = D2(i, j); //Be
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			af = vectorToAAFR(oB(i, j));
			Bf(i, j) = af;
		}
	}
	for (int i = 0; i < n; i++) {
		b_i_ob(i, obn, bi);
		af = vectorToAAFR(np, bi);
		bf[i] = af;
	}
	aafrvector wf = Bf * bf;
	for (int i = 0; i < n; i++) {
		w[i] = wf[i].reduce();
	}
	t = difftime(clock(), start) / CLOCKS_PER_SEC;
}

//%============================================================================
//%
//% Degrauwe's method for general parametric interval linear systems
//%
//%============================================================================

bool
//%----------------------------------------------------------------------------
//% Main to verify Degrauwe method for solving problems. In this approach 
//% A_e=0, A_e normally corresponds to the last element in reduced affine forms
//%----------------------------------------------------------------------------
Degrauwe1Order(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w,
	real& t)
{
	AAFR af;
	dmatrix	R, Bi, Ci, C2, C3, D2;
	aafrmatrix Bf;
	aafrvector bf;

	int n = b1[0].size();
	int K = p1.size();
	dmatrix Z(n, n), I(n, n);
	ivector p(K);
	omatrix oB(n, n);
	omvector M(K + 1), B(K + 1);
	ovector b(K + 1), bb(K + 1), obn(n);

	ktransform(M1, b1, p1, M, b, p); // the system is transformed so that all [p]=[-1,1]

	R = M[0];
	if (inv(R)) {
		dvector el = dvector(K + 2);
		dvector el2 = dvector(K + 1);
		el.fill_in(0.0);
		oB.fill_in(el);
		obn.fill_in(el2);
		Ci = dmatrix(n, n);
		w = ivector(n);
		Z.fill_in(0.0);
		Bf = aafrmatrix(n, n);
		bf = aafrvector(n);
		//% The estimation of the inverse has the form: B0+sum_{i=1}^KB_i*e_i+B_e*e_e
		//% B0=I, B_i=-A_i, B_e=(Z^2*(I-Z)^{-1})
		Unit(I);
		B[0] = I; // B0 = I
		for (int i = 0; i < n; i++) {
			bb[0] = R * b[0];
			obn[i][0] = bb[0][i];
			for (int j = 0; j < n; j++) {
				oB(i, j)[0] = B[0](i, j);
			}
		}
		for (int k = 1; k < K + 1; k++) {
			B[k] = (R * M[k]) * (-1.0); //% Bk = -(R * Ai)
			bb[k] = R * b[k];
			for (int i = 0; i < n; i++) {
				obn[i][k] = bb[k][i];
				for (int j = 0; j < n; j++) {
					oB(i, j)[k] = B[k](i, j);
				}
			}
			mag(B[k], Ci);
			Z = Z + Ci; //% Z=sum_{i=1}^np |R*Ai|
		}

		C2 = Z * Z; //% C^2
		C3 = I - Z; //% C3=I-C
		if (inv(C3)) { //% (I-C)^{-1}
			D2 = C2 * C3; //% B_e=Z^2*(I-Z)^{-1} - part enclosing the remainder of the Neumann series
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					oB(i, j)[K + 1] = D2(i, j); //Be
				}
			}
			//% transforming to affine forms
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					af = vectorToAAFR(oB(i, j));
					Bf(i, j) = af;
				}
			}
			for (int i = 0; i < n; i++) {
				//% if the vector of parameters in the matrix
				af = vectorToAAFR(1, obn[i]);
				bf[i] = af;
			}
			aafrvector wf = Bf * bf;
			for (int i = 0; i < n; i++) {
				w[i] = wf[i].reduce();
			}
			return true;
		}
		else {
			return false;
		}
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Main to verify 2-order Degrauwe method
//%----------------------------------------------------------------------------
Degrauwe2Order(const omvector& M1,
	const ovector& b1,
	const ivector& p1,
	ivector& w,
	real& t)
{
	int n = b1[0].size();
	int K = p1.size();
	aafrmatrix Bf(n, n);
	aafrvector bf(n);
	dvector	el, el2;
	dmatrix A2(n, n), A2i(n, n), A3(n, n), A4(n, n), Zi(n, n);
	dmatrix	Z(n, n), I(n, n), R(n, n), C2(n, n), C3(n, n), D2(n, n), D3(n, n), Fi(n, n), Fii(n, n);
	ivector p(K);
	omatrix oB(n, n);
	omvector M(K + 1), B(K + 2);
	ovector b(K + 1), bb(K + 1), obn(n);

	ktransform(M1, b1, p1, M, b, p); //% the system is transformed so that all [p]=[-1,1]
	R = M[0];
	if (inv(R)) {
		w = ivector(n);
		dvector el = dvector(K + 2);
		dvector el2 = dvector(K + 1);
		el.fill_in(0.0);
		oB.fill_in(el);
		obn.fill_in(el2);
		//% left pre-conditioning
		for (int k = 0; k < K + 1; k++) {
			M[k] = R * M[k];
			bb[k] = R * b[k];
		}
		Z.fill_in(0.0);
		A2.fill_in(0.0);
		A3.fill_in(0.0);
		//% The estimation of the inverse has the form: B0+sum_{i=1}^KB_i*e_i+B_e*e_e
		//% B0=I+0.5A^2_i, B_i=-A_i, B_e=(Z^3*(I-Z)^{-1})
		Unit(I);
		for (int k = 1; k < K + 1; k++) {
			mag(M[k], D3); // |Ai|
			Z = Z + D3; // sum(|Ai|)
			D2 = M[k] * M[k]; // Ai^2
			A2 = A2 + D2; // sum(Ai^2)
			mag(D2, D3); // |Ai^2|
			A3 = A3 + D3; // sum(|Ai^2|
		}
		B[0] = I + A2 * 0.5; //% B0 = I
		for (int k = 1; k < K + 1; k++) {
			B[k] = M[k] * (-1.0); //% Bk = -Bi
		}
		A4.fill_in(0.0);
		for (int k = 1; k < K + 1; k++) {
			for (int kk = k + 1; kk < K + 1; kk++) {
				D2 = M[k] * M[kk] + M[kk] * M[k];
				mag(D2, D3);
				A4 = A4 + D3; // sum(|Ai*Aj+Aj*Ai|)
			}
		}
		C3 = I - Z;
		inv(C3);
		B[K + 1] = A3 * 0.5 + A4 + Z * Z * Z * C3;
		for (int k = 0; k < K + 1; k++) {
			for (int i = 0; i < n; i++) {
				obn[i][k] = bb[k][i];
				for (int j = 0; j < n; j++) {
					oB(i, j)[k] = B[k](i, j);
				}
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				oB(i, j)[K + 1] = B[K + 1](i, j); //Be
			}
		}
		//% brakuje uwzglednienia Be, TRZEBA KONIECZNIE DODAC!!!!!
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Bf(i, j) = vectorToAAFR(oB(i, j));
			}
		}
		for (int i = 0; i < n; i++) {
			//% if the vector of parameters in the matrix
			bf[i] = vectorToAAFR(1, obn[i]);
		}
		aafrvector wf = Bf * bf;
		for (int i = 0; i < n; i++) {
			w[i] = wf[i].reduce();
		}
		return true;
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Modified Expansion Method with residual correction
//% with additional modification of Milan
//% This looks like not finished!!!
//%----------------------------------------------------------------------------
DegrauveMM(int m, //% order of the method
	const aafrmatrix& A, //% revised affine matrix
	const aafrvector& b, //% revised affine right-hand side vector
	aafrvector& w)		 //% revised affine solution vector
{
	int n = A.num_rows();
	aafrmatrix AA(n, n), V(n, n);
	aafrvector v(n), d(n);
	dmatrix R(n, n), I(n, n);
	dvector bm(n), x0(n);

	mid(A, R); //% midpoint matrix
	mid(b, bm); //% midpoint vector
	if (inv(R)) { //% midpoint inverse
		x0 = R * bm; //% midpoint solution
		AA = R * A; //% left pre-conditioning
		v = R * (b - A * x0); //% left pre-conditioning and residual correction
		Unit(I);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				V(i, j) = I(i, j) - AA(i, j);
			}
		}
		d = v;
		for (int i = 0; i < m; i++) {
			d = V * d;
		}
		return true;
	}
	return false;
}