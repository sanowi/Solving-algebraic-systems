#include "../utils/stdafx.h"
#include <conio.h>
#include "../vector_matrix/vector_matrix_base.h"
#include "../affine/aafr.h"
#include "../utils/inverse.h"
#include "../utils/randnum.h"
#include "../utils/gausselim.h"
#include "imllib.h"

//http://math.nist.gov/iml++/qmr.h.txt

int
CG(const aafmatrix &A, ivector &w, const aafvector &b, dmatrix &M, int &max_iter, real &tol)
{
	int n = b.size();
	real resid;
	AAF alpha, beta, rho, rho_1;
	dmatrix R(n, n);
	aafvector p, z, q, x(n), bb(n);
	aafmatrix AA(n, n);
	
	AAF::setlast(0);
	x.fill_in(0.0);
	w = ivector(n);

	mid(A, R);
	if (inv(R)) {
		AA = R*A;
		bb = R*b;

		aafvector r = bb - AA * x; // initial

		resid = 0.0;
		for (int i = 0; i < n; i++) {
			resid += r[i].reduce().rad();
		}
		std::cout << "resid = " << resid << std::endl;
		if (resid <= tol) {
			for (int i = 0; i < n; i++) {
				w[i] = x[i].reduce();
			}
			tol = resid;
			max_iter = 0;
			return 0;
		}

		Unit(M);

		for (int i = 1; i <= max_iter; i++) {
			z = r;
			AAF s = 0.0;
			for (int j = 0; j < n; j++) { // dot(r, z);
				s = s + r[j] * z[j];
			}
			rho = s;

			if (i == 1) {
				p = z;
			}
			else {
				beta = rho / rho_1;
				p = z + beta * p;
			}

			q = AA * p;

			s = 0.0;
			for (int j = 0; j < n; j++) {
				s = s + p[j] * q[j];
			}
			alpha = rho / s; // dot(p, q);

			x = x + alpha * p;
			r = r - alpha * q;

			resid = 0.0;
			for (int i = 0; i < n; i++) {
				resid += r[i].reduce().rad();
			}

			if (resid <= tol) {
				for (int i = 0; i < n; i++) {
					w[i] = x[i].reduce();
				}
				tol = resid;
				max_iter = 0;
				return 0;
			}

			rho_1 = rho;
		}

		tol = resid;
		return 1;
	}
	return 1;
}

int 
CG(const dmatrix &A, dvector &x, const dvector &b,dmatrix &M, int &max_iter, real &tol)
{
	int n = x.size();
	real resid;
	dvector p, z, q;
	dvector alpha(1), beta(1), rho(1), rho_1(1);

	real normb = norm(b);
	dvector r = b - A*x; // initial

	if (normb == 0.0) {
		normb = 1.0;
	}

	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	Unit(M);

	for (int i = 1; i <= max_iter; i++) {
		//z = M.solve(r);
		gw_gausse(M, r, z);
		z = r;
		real s = 0.0;
		for (int j = 0; j < n; j++){ // dot(r, z);
			s = s + r[j] * z[j];  
		}
		rho[0] = s;

		if (i == 1) {
			p = z;
		}
		else {
			beta[0] = rho[0] / rho_1[0];
			p = z + beta[0] * p;
		}

		q = A * p;

		s = 0.0;
		for (int j = 0; j < n; j++) {
			s = s + p[j] * q[j];
		}
		alpha[0] = rho[0] / s; // dot(p, q);

		x += alpha[0] * p;
		r -= alpha[0] * q;

		if ((resid = norm(r) / normb) <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		rho_1[0] = rho[0];
	}

	tol = resid;
	return 1;
}

//int
//BiCG(const Matrix &A, Vector &x, const Vector &b,
//const Preconditioner &M, int &max_iter, Real &tol)
//{
//	Real resid;
//	Vector rho_1(1), rho_2(1), alpha(1), beta(1);
//	Vector z, ztilde, p, ptilde, q, qtilde;
//
//	Real normb = norm(b);
//	Vector r = b - A * x;
//	Vector rtilde = r;
//
//	if (normb == 0.0)
//		normb = 1;
//
//	if ((resid = norm(r) / normb) <= tol) {
//		tol = resid;
//		max_iter = 0;
//		return 0;
//	}
//
//	for (int i = 1; i <= max_iter; i++) {
//		z = M.solve(r);
//		ztilde = M.trans_solve(rtilde);
//		rho_1(0) = dot(z, rtilde);
//		if (rho_1(0) == 0) {
//			tol = norm(r) / normb;
//			max_iter = i;
//			return 2;
//		}
//		if (i == 1) {
//			p = z;
//			ptilde = ztilde;
//		}
//		else {
//			beta(0) = rho_1(0) / rho_2(0);
//			p = z + beta(0) * p;
//			ptilde = ztilde + beta(0) * ptilde;
//		}
//		q = A * p;
//		qtilde = A.trans_mult(ptilde);
//		alpha(0) = rho_1(0) / dot(ptilde, q);
//		x += alpha(0) * p;
//		r -= alpha(0) * q;
//		rtilde -= alpha(0) * qtilde;
//
//		rho_2(0) = rho_1(0);
//		if ((resid = norm(r) / normb) < tol) {
//			tol = resid;
//			max_iter = i;
//			return 0;
//		}
//	}
//
//	tol = resid;
//	return 1;
//}

//int
//QMR(const Matrix &A, Vector &x, const Vector &b, const Preconditioner1 &M1,
//const Preconditioner2 &M2, int &max_iter, Real &tol)
//{
//	Real resid;
//
//	Vector rho(1), rho_1(1), xi(1), gamma(1), gamma_1(1), theta(1), theta_1(1);
//	Vector eta(1), delta(1), ep(1), beta(1);
//
//	Vector r, v_tld, y, w_tld, z;
//	Vector v, w, y_tld, z_tld;
//	Vector p, q, p_tld, d, s;
//
//	Real normb = norm(b);
//
//	r = b - A * x;
//
//	if (normb == 0.0)
//		normb = 1;
//
//	if ((resid = norm(r) / normb) <= tol) {
//		tol = resid;
//		max_iter = 0;
//		return 0;
//	}
//
//	v_tld = r;
//	y = M1.solve(v_tld);
//	rho(0) = norm(y);
//
//	w_tld = r;
//	z = M2.trans_solve(w_tld);
//	xi(0) = norm(z);
//
//	gamma(0) = 1.0;
//	eta(0) = -1.0;
//	theta(0) = 0.0;
//
//	for (int i = 1; i <= max_iter; i++) {
//
//		if (rho(0) == 0.0)
//			return 2;                        // return on breakdown
//
//		if (xi(0) == 0.0)
//			return 7;                        // return on breakdown
//
//		v = (1. / rho(0)) * v_tld;
//		y = (1. / rho(0)) * y;
//
//		w = (1. / xi(0)) * w_tld;
//		z = (1. / xi(0)) * z;
//
//		delta(0) = dot(z, y);
//		if (delta(0) == 0.0)
//			return 5;                        // return on breakdown
//
//		y_tld = M2.solve(y);               // apply preconditioners
//		z_tld = M1.trans_solve(z);
//
//		if (i > 1) {
//			p = y_tld - (xi(0) * delta(0) / ep(0)) * p;
//			q = z_tld - (rho(0) * delta(0) / ep(0)) * q;
//		}
//		else {
//			p = y_tld;
//			q = z_tld;
//		}
//
//		p_tld = A * p;
//		ep(0) = dot(q, p_tld);
//		if (ep(0) == 0.0)
//			return 6;                        // return on breakdown
//
//		beta(0) = ep(0) / delta(0);
//		if (beta(0) == 0.0)
//			return 3;                        // return on breakdown
//
//		v_tld = p_tld - beta(0) * v;
//		y = M1.solve(v_tld);
//
//		rho_1(0) = rho(0);
//		rho(0) = norm(y);
//		w_tld = A.trans_mult(q) - beta(0) * w;
//		z = M2.trans_solve(w_tld);
//
//		xi(0) = norm(z);
//
//		gamma_1(0) = gamma(0);
//		theta_1(0) = theta(0);
//
//		theta(0) = rho(0) / (gamma_1(0) * beta(0));
//		gamma(0) = 1.0 / sqrt(1.0 + theta(0) * theta(0));
//
//		if (gamma(0) == 0.0)
//			return 4;                        // return on breakdown
//
//		eta(0) = -eta(0) * rho_1(0) * gamma(0) * gamma(0) /
//			(beta(0) * gamma_1(0) * gamma_1(0));
//
//		if (i > 1) {
//			d = eta(0) * p + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * d;
//			s = eta(0) * p_tld + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * s;
//		}
//		else {
//			d = eta(0) * p;
//			s = eta(0) * p_tld;
//		}
//
//		x += d;                            // update approximation vector
//		r -= s;                            // compute residual
//
//		if ((resid = norm(r) / normb) <= tol) {
//			tol = resid;
//			max_iter = i;
//			return 0;
//		}
//	}
//
//	tol = resid;
//	return 1;                            // no convergence
//}



//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************


//template < class Matrix, class Vector >
//void
//Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[])
//{
//	Vector y(s);
//
//	// Backsolve:  
//	for (int i = k; i >= 0; i--) {
//		y(i) /= h(i, i);
//		for (int j = i - 1; j >= 0; j--)
//			y(j) -= h(j, i) * y(i);
//	}
//
//	for (int j = 0; j <= k; j++)
//		x += v[j] * y(j);
//}
//
//
//template < class Real >
//Real
//abs(Real x)
//{
//	return (x > 0 ? x : -x);
//}
//
//
//template < class Operator, class Vector, class Preconditioner,
//class Matrix, class Real >
//	int
//	GMRES(const Operator &A, Vector &x, const Vector &b,
//	const Preconditioner &M, Matrix &H, int &m, int &max_iter,
//	Real &tol)
//{
//	Real resid;
//	int i, j = 1, k;
//	Vector s(m + 1), cs(m + 1), sn(m + 1), w;
//
//	Real normb = norm(M.solve(b));
//	Vector r = M.solve(b - A * x);
//	Real beta = norm(r);
//
//	if (normb == 0.0)
//		normb = 1;
//
//	if ((resid = norm(r) / normb) <= tol) {
//		tol = resid;
//		max_iter = 0;
//		return 0;
//	}
//
//	Vector *v = new Vector[m + 1];
//
//	while (j <= max_iter) {
//		v[0] = r * (1.0 / beta);    // ??? r / beta
//		s = 0.0;
//		s(0) = beta;
//
//		for (i = 0; i < m && j <= max_iter; i++, j++) {
//			w = M.solve(A * v[i]);
//			for (k = 0; k <= i; k++) {
//				H(k, i) = dot(w, v[k]);
//				w -= H(k, i) * v[k];
//			}
//			H(i + 1, i) = norm(w);
//			v[i + 1] = w * (1.0 / H(i + 1, i)); // ??? w / H(i+1, i)
//
//			for (k = 0; k < i; k++)
//				ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));
//
//			GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
//			ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
//			ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));
//
//			if ((resid = abs(s(i + 1)) / normb) < tol) {
//				Update(x, i, H, s, v);
//				tol = resid;
//				max_iter = j;
//				delete[] v;
//				return 0;
//			}
//		}
//		Update(x, i - 1, H, s, v);
//		r = M.solve(b - A * x);
//		beta = norm(r);
//		if ((resid = beta / normb) < tol) {
//			tol = resid;
//			max_iter = j;
//			delete[] v;
//			return 0;
//		}
//	}
//
//	tol = resid;
//	delete[] v;
//	return 1;
//}
//
//
//#include <math.h> 
//
//
//template<class Real>
//void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
//{
//	if (dy == 0.0) {
//		cs = 1.0;
//		sn = 0.0;
//	}
//	else if (abs(dy) > abs(dx)) {
//		Real temp = dx / dy;
//		sn = 1.0 / sqrt(1.0 + temp*temp);
//		cs = temp * sn;
//	}
//	else {
//		Real temp = dy / dx;
//		cs = 1.0 / sqrt(1.0 + temp*temp);
//		sn = temp * cs;
//	}
//}
//
//
//template<class Real>
//void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
//{
//	Real temp = cs * dx + sn * dy;
//	dy = -sn * dx + cs * dy;
//	dx = temp;
//}