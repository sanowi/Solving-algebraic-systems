#include "../vector_matrix/vector_matrix_base.h"
#include "../utils/inverse.h"
#include "../utils/gausselim.h"
#include "../affine/aafr.h"

bool
//%----------------------------------------------------------------------------
//% Method for solving systems of the form: A^T * D * A * x = b,
//% where A is an (n x n) real matrix,
//% (Caution: A is an argument not A^T)
//% the method involves affine arithmetic
//% (Caution!!!: It works only for square A, which means that
//% only statically determined trusses can be solved using this method)
//%----------------------------------------------------------------------------
mneums2(const dmatrix& A,
	const ivector& p,
	const ivector& q,
	ivector& w)
{
	int	np = p.size(), nq = q.size(); // q to wlasciwie to samo co b
	aafrvector pv(np), qv(nq), x(nq), y(np), z(np), aw(np);
	cvector ind(np), indx(nq);
	dmatrix A0 = A, R, S, Q, R2; // , AT = A.transpose()

	w = ivector(nq);
	for (int i = 0; i < np; i++) {
		pv[i] = AAFR(p[i]);
		aw[i] = 0.0;
	}
	for (int i = 0; i < nq; i++) {
		qv[i] = AAFR(q[i]);
	}
	R = A;
	gw_gausse(A, qv, aw);
	for (int i = 0; i < np; i++) {
		y[i] = aw[i] / pv[i];
	}
	R = A.transpose();
	gw_gausse(R, y, x);
	for (int i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}

bool
//%----------------------------------------------------------------------------
//% Method for solving systems of the form: A^T * D * A * x = b,
//% where A is an (n x n) real matrix,
//% (Caution: A is an argument not A^T)
//% the method involves affine arithmetic
//% here: Ri = inv(A), Si = inv(A^T)
//% (Caution: In case of mechanical problems: it works only for square A,
//% which means that only statically determine trusses can be solved using this method)
//%----------------------------------------------------------------------------
mneumss(const dmatrix& A,
	const ivector& p,
	const ivector& q,
	ivector& w)
{
	int i, j, np = p.size(), nq = q.size(); //% q is acutally the same as b
	AAFR aw;
	aafrvector pv(np), qv(nq), y(np), x(nq), z(np);
	cvector idx(2), idx1(1), ind(np);
	dvector cfs(2), cfs1(1);
	cvector indx(nq);
	dmatrix A0 = A;
	dmatrix	AT = A.transpose(), R, S, Q, R2;
	imatrix Ri, Si;
	aafrmatrix Rf, Sf;
	
	w = ivector(nq);
	for (i = 0; i < np; i++) {
		pv[i] = AAFR(p[i]);
	}
	for (i = 0; i < nq; i++) {
		qv[i] = AAFR(q[i]);
	}
	R = A;
	Ri = imatrix(R.num_rows(), R.num_cols());
	to_int(R, Ri);
	inv(Ri); // wyglada na to, ze odwrocenie macierzy jest najkorzystniejsze
	Rf = aafrmatrix(Ri.num_rows(), Ri.num_cols());
	Sf = aafrmatrix(Ri.num_cols(), Ri.num_rows());
	S = A.transpose();
	Si = imatrix(R.num_cols(), R.num_rows());
	to_int(S, Si);
	inv(Si);
	for (int i = 0; i < R.num_rows(); i++) {
		for (j = 0; j < R.num_cols(); j++) {
			Rf(i, j) = AAFR(Ri(i, j));
			Sf(j, i) = AAFR(Si(j, i));
		}
	}
	for (i = 0; i < np; i++) { //% tu moze nie od i = 0, tylko od i = nq + 1, czy jakos tak
		//% bo rozwiazanie stanowi koncowka wektora, a niepoczatek, porownac z artykulem Popovej
		aw = 0;
		for (j = 0; j < nq; j++) {
			if (!(Sf(i, j) == 0) && !(qv[j] == 0)) {
				aw = aw + Sf(i, j) * qv[j]; // aw - rozwiazanie ukladu A^T * aw = qv - uklad niedookreslony
			}
		}
		y[i] = aw / pv[i]; //* (1.0/pv[i]); // y - rozwiazanie ukladu D * y = aw - uklad kwadratowy
	}
	for (i = 0; i < nq; i++) {
		aw = 0;
		for (j = 0; j < np; j++) {
			if (!(Rf(i, j) == 0) && !(y[j] == 0)) {
				aw = aw + Rf(i, j) * y[j]; // aw - rozwiazanie ukladu A * aw = y - uklad nadokreslony
			}
		}
		x[i] = aw;
	}
	for (i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}

bool
//%----------------------------------------------------------------------------
//% Method for solving systems of the form: A^T * D * A * x = b,
//% where A is an (n x n) real matrix,
//% (Caution: A is an argument not A^T)
//% the method involves affine arithmetic
//% here: Ri = inv(A), Si = Ri^T
//% (Caution: In case of mechanical problems: it works only for square A,
//% which means that only statically determine trusses can be solved using this method)
//%----------------------------------------------------------------------------
mneumsss(const dmatrix& A,
	const ivector& p,
	const ivector& q,
	ivector& w)
{
	int np = p.size(), nq = q.size(); //% q is actually the same as b
	aafrvector pv(np), qv(nq), y(np), x(nq), z(np); 
	AAFR aw;
	cvector idx(2), idx1(1), ind(np);
	dvector cfs(2), cfs1(1);
	cvector indx(nq);
	dmatrix A0 = A, AT = A.transpose(), R, S, Q, R2;
	imatrix Ri, Si;
	aafrmatrix Rf, Sf;
	
	w = ivector(nq);
	for (int i = 0; i < np; i++) {
		pv[i] = AAFR(p[i]);
	}
	for (int i = 0; i < nq; i++) {
		qv[i] = AAFR(q[i]);
	}
	R = A;
	Ri = imatrix(R.num_rows(), R.num_cols());
	to_int(R, Ri);
	inv(Ri); // wyglada na to, ze odwrocenie macierzy jest najkorzystniejsze
	Si = Ri.transpose();
	Rf = aafrmatrix(Ri.num_rows(), Ri.num_cols());
	Sf = aafrmatrix(Si.num_rows(), Si.num_cols());
	for (int i = 0; i < Rf.num_rows(); i++) {
		for (int j = 0; j < Rf.num_cols(); j++) {
			Rf(i, j) = AAFR(Ri(i, j));
			Sf(j, i) = AAFR(Si(j, i));
		}
	}
	for (int i = 0; i < np; i++) { //% tu moze nie od i = 0, tylko od i = nq + 1, czy jakos tak
		// bo rozwiazanie stanowi koncowka wektora, a niepoczatek, porownac z artykulem Popovej
		aw = 0;
		for (int j = 0; j < nq; j++) {
			if (!(Sf(i, j) == 0) && !(qv[j] == 0)) {
				aw = aw + Sf(i, j) * qv[j]; // aw - rozwiazanie ukladu A^T * aw = qv - uklad niedookreslony
			}
		}
		y[i] = aw * (1.0 / pv[i]); // y - rozwiazanie ukladu D * y = aw - uklad kwadratowy
	}
	for (int i = 0; i < nq; i++) {
		aw = 0;
		for (int j = 0; j < np; j++) {
			if (!(Rf(i, j) == 0) && !(y[j] == 0)) {
				aw = aw + Rf(i, j) * y[j]; // aw - rozwiazanie ukladu A * aw = y - uklad nadokreslony
			}
		}
		x[i] = aw;
	}
	for (int i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}

bool
//%----------------------------------------------------------------------------
//% solving systems of the form: A^T * D * A * x = b
//% A	- real matrix
//% this method solves the following systems:
//% (A^T) * z = b (<=> z = (A^T)^-1 * b)
//% D * y = z (<=> y = D^{-1}z
//% A * x = y
//% Remark: A can be rectangular!!!
//%----------------------------------------------------------------------------
mneumr(const dmatrix& A,
	const ivector& p,
	const ivector& q,
	ivector& w)
{
	//% q = b, np liczba wierszy, nq - liczba kolumn, dla niekwadratowych np > nq
	int	m = A.num_rows(), n = A.num_cols(), np = p.size(), nq = q.size();
	AAFR aw, s0;
	aafrvector	pv(np), qv(nq), y(np), y2(np), x(nq), b1(m + n); //% b1 - extended right hand vector
	dmatrix	B(m + n, m + n), Id, G, R;
	
	for (int i = 0; i < np; i++) { //% left hand interval vector
		pv[i] = AAFR(p[i]);
	}
	for (int i = 0; i < nq; i++) { //% right hand interval vector
		qv[i] = AAFR(q[i]);
	}
	w = ivector(nq);
	B.fill_in(0);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			B(i, j) = A(i, j);
			B(j + m, i + n) = A(i, j);
		}
	}
	for (int i = 0; i < m; i++) {
		B(i, i + n) = -1;
	}
	R = B; //% B = [A^T -I]
	//%     [0    A]
	inv(R); //% R = B^{-1}
	for (int i = 0; i < np; i++) {
		aw = 0;
		for (int j = 0; j < nq; j++) { //% solving underdetermined system
			if (R(i + n, j + m) != 0 && !(qv[j] == 0)) { // bylo R
				aw = aw + R(i + n, j + m) * qv[j]; // bylo R
			}
		}
		y[i] = aw / pv[i]; //% (1.0/pv[i]);
	}
	//% solution of the overdermined system
	for (int i = 0; i < nq; i++) { //% takes lot of time
		aw = 0;
		for (int j = 0; j < np; j++) {
			if (R(i, j) != 0 && !(y[j] == 0)) {
				aw = aw + R(i, j) * y[j];
			}
		}
		x[i] = aw;
	}
	for (int i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}

bool
//%----------------------------------------------------------------------------
//% Solving systems of the form: A^T * D * A * x = b
//% A	- real matrix
//% this method solves the following systems:
//% (A^T) * z = b (<=> z = (A^T)^-1 * b)
//% D * y = z (<=> y = D^{-1}z
//% A * x = y
//% Remark: A can be rectangular!!!
//%----------------------------------------------------------------------------
mneumrr(const dmatrix& A, 
	const ivector& p, 
	const ivector& q, 
	ivector& w)
{
	//% q = b, np number of rows, nq - number of columns, for rectangular systems np > nq
	int	m = A.num_rows(), n = A.num_cols(), np = p.size(), nq = q.size();
	AAFR aw, s0;
	aafrvector	pv(np), qv(nq), y(np), y2(np), x(nq), b1(m + n); //% b1 - extended right hand vector
	dmatrix	B(m + n, m + n), Id, G, R;
	imatrix Ri(m + n, m + n);
	aafrmatrix	Rf;

	for (int i = 0; i < np; i++) { //% left hand interval vector
		pv[i] = AAFR(p[i]);
	}
	for (int i = 0; i < nq; i++) { //% right hand interval vector
		qv[i] = AAFR(q[i]);
	}
	w = ivector(nq);
	B.fill_in(0);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			B(i, j) = A(i, j);
			B(j + m, i + n) = A(i, j);
		}
	}
	for (int i = 0; i < m; i++) {
		B(i, i + n) = -1;
	}
	R = B; //% B = [A^T -I]
	       //%     [0    A]
	to_int(R, Ri);
	inv(Ri); // R = B^{-1}
	Rf = aafrmatrix(Ri.num_rows(), Ri.num_cols());
	for (int i = 0; i < R.num_rows(); i++) {
		for (int j = 0; j < R.num_cols(); j++) {
			Rf(i, j) = AAFR(Ri(i, j));
		}
	}
	for (int i = 0; i < np; i++) {
		aw = 0;
		for (int j = 0; j < nq; j++)  { //% solving underdetermined system
			if (!(Rf(i + n, j + m) == 0) && !(qv[j] == 0)) {
				aw = aw + Rf(i + n, j + m) * qv[j];
			}
		}
		y[i] = aw  * (1.0 / pv[i]);
	}
	//% solution of the overestimated system (this generates the error 
	//% since we actally obtain a solution which minimizes quadratic error
	for (int i = 0; i < nq; i++) {
		aw = 0;
		for (int j = 0; j < np; j++) {
			if (!(Rf(i, j) == 0) && !(y[j] == 0)) {
				aw = aw + Rf(i, j) * y[j];
			}
		}
		x[i] = aw;
	}
	for (int i = 0; i < nq; i++) {
		w[i] = x[i].reduce();
	}
	return true;
}