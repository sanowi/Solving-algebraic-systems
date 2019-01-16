#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../affine/aafr.h"

const real NEUMHL_EPS = 1.0e-10;

bool
//%----------------------------------------------------------------------------
//% Iterative method for solving systems of the form
//% K + B*D(p)*A = a + F*b,
//% where B and A are real square matrices (!!!)
//% Hladik's modification of the method
//% from paper "Linear systems with large uncertainties..."
//% by Neumaier and Pownuk
//%----------------------------------------------------------------------------
neumhl(const dmatrix& K, const dmatrix& A, const dmatrix& B, const imatrix& D, const dmatrix& F, const dvector& a, const ivector& b, ivector& w)
{
	int	n = A.num_rows();; // ne - number of edges, dof = ne + 1
	dmatrix BB = B;
	dmatrix AA = A;
	if (inv(BB) && inv(AA)) {
		dmatrix KK = BB * K * AA;
		ivector bb = BB * (a + F * b);
		imatrix CC = KK + D;
		inv(CC);
		w = AA* CC * bb;
		return true;
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Method for solving systems of the form: B * D * A * x = b,
//% where A is an (n x n) real matrix,
//% (Caution: A is an argument not A^T)
//% the method involves affine arithmetic
//% (Caution!!!: It works only for square A, which means that
//% only statically determine trusses can be solved using this method)
//%----------------------------------------------------------------------------
neums(const dmatrix& B, const dmatrix& A, const ivector& p, const ivector& q, ivector& w)
{
	int	np = p.size(), nq = q.size();
	AAFR aw;
	aafrvector	pv(np), qv(nq), x(nq), y(np), z(np);
	cvector	ind(np), indx(nq);
	dmatrix	A0 = A, R, S, R2;

	w = ivector(nq);
	for (int i = 0; i < np; i++) {
		pv[i] = AAFR(p[i]);
	}
	for (int i = 0; i < nq; i++) {
		qv[i] = AAFR(q[i]);
	}
	R = B;
	S = A;
	if (inv(R) && inv(S)) { 
		for (int i = 0; i < np; i++) {
			aw = 0.0;
			for (int j = 0; j < nq; j++) {
				if (S(i, j) != 0.0 && !(qv[j] == 0.0)) {
					aw = aw + S(i, j) * qv[j]; // 
				}
			}
			y[i] = aw / pv[i]; //(1.0/pv[i]); // y - rozwiazanie ukladu D * y = aw - uklad kwadratowy
		}
		for (int i = 0; i < nq; i++) {
			aw = 0.0;
			for (int j = 0; j < np; j++) {
				if (R(i, j) != 0.0 && !(y[j] == 0.0)) {
					aw = aw + R(i, j) * y[j]; // 
				}
			}
			x[i] = aw;
		}
		for (int i = 0; i < nq; i++) {
			w[i] = x[i].reduce();
		}
		return true;
	}
	return false;
}

bool
//%----------------------------------------------------------------------------
//% Method for solving systems of the form: A^T * D * A * x = b,
//% where A is an (n x n) real matrix,
//% (Caution: A is an argument not A^T)
//% the method involves affine arithmetic
//% (Caution!!!: It works only for square A, which means that
//% only statically determine trusses can be solved using this method)
//%----------------------------------------------------------------------------
neums(const dmatrix& A, const ivector& p, const ivector& q, ivector& w)
{
	int	np = p.size(), nq = q.size();
	AAFR aw;
	aafrvector	pv(np), qv(nq), x(nq), y(np), z(np);
	cvector	ind(np), indx(nq);
	dmatrix	A0 = A, R, S, R2;

	w = ivector(nq);
	for (int i = 0; i < np; i++) {
		pv[i] = AAFR(p[i]);
	}
	for (int i = 0; i < nq; i++) {
		qv[i] = AAFR(q[i]);
	}
	R = A;
	//S = AT;
	inv(R); // wyglada na to, ze odwrocenie macierzy jest najkorzystniejsze
	//inv( S );
	S = R.transpose();
	for (int i = 0; i < np; i++) { // tu moze nie od i = 0, tylko od i = nq + 1, czy jakos tak
		// bo rozwiazanie stanowi koncowka wektora, a niepoczatek, porownac z artykulem Popovej
		aw = 0.0;
		for (int j = 0; j < nq; j++) {
			if (S(i, j) != 0.0 && !(qv[j] == 0.0)) {
				aw = aw + S(i, j) * qv[j]; // aw - rozwiazanie ukladu A^T * aw = qv - uklad niedookreslony
			}
		}
		y[i] = aw / pv[i]; //(1.0/pv[i]); // y - rozwiazanie ukladu D * y = aw - uklad kwadratowy
	}
	for (int i = 0; i < nq; i++) {
		aw = 0.0;
		for (int j = 0; j < np; j++) {
			if (R(i, j) != 0.0 && !(y[j] == 0.0)) {
				aw = aw + R(i, j) * y[j]; // aw - rozwiazanie ukladu A * aw = y - uklad nadokreslony
			}
		}
		x[i] = aw;
	}
	for (int i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}

