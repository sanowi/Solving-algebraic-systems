#include <iostream>
#include <math.h>
#include <time.h>
#include "config.h"
#include "inverse.h"
#include "mism.h"
#include "TrussStructure.h"
#include "aafr.h"
#include "rohn.h"

bool Rohn(	const ivector& p, 
			const omatrix& oA, 
			const ovector& ob, 
			ivector& w)
{
	int			n =	ob.size(), np = p.size();
	REAL		eps = 1.0e-16;
	dvector		x0(n), d(n), pm(np), bm(n), d1(n), pom(n), dpom(n);
	dmatrix		R(n, n), M(n, n), I(n, n), Gm(n, n);
	ivector		z(n), x(n), x1(n);
    imatrix		C(n, n), Rm(n, n);
	//-----------------------
	mid(p, pm);
	bm = v_from_dep(pm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	if (inv(R)) 
	{
		x0 = R * bm; // to jest x*
		c_matrix(R, p, oA, C); // C = mid(A(p))^{-1}*A(p)
		z_vector(x0, R, p, oA, ob, z); // g = b = mid(A(p)^{-1}*b(q)
		Unit(I);
		Rm = I - C; // G = Rm
		mag(Rm, Gm);
		d1.fill_in(0);
		x1.fill_in(0);
		do 
		{
			x = x1;
			d = d1;
			x1 = Rm * x + z;
			mag(x1 - x, pom);
			d1 = Gm * d + pom + eps;
			mag(d1 - d, dpom);
		}
		while (!(dpom < eps));
		w = x0 + d1 * interval(-1.0, 1.0);
		return true;
	}
	return false;
}

bool Rohn(	const ivector& p, 
			const ivector& q, 
			const omatrix& oA, 
			const ovector& ob, 
			ivector& w)
{
	int			n =	ob.size(), np = p.size(), nq = q.size();
	REAL		eps = 1.0e-16;
	dvector		x0(n), d(n), pm(np), qm(nq), bm(n), d1(n), pom(n), dpom(n);
	dmatrix		R(n, n), M(n, n), I(n, n), Gm(n, n);
	ivector		x(n), x1(n), z(n);
    imatrix		C(n, n), Rm(n, n);
	//-------------------
	mid(p, pm);
	mid(q, qm);
	bm = v_from_dep(qm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	if (inv(R)) 
	{
		x0 = R * bm; // to jest x*
		c_matrix(R, p, oA, C); // C = R*A(p)
		z_vector(x0, R, p, q, oA, ob, z); // g = b = R*(b - A*x0)
		Unit(I);
		Rm = I - C; // G = Rm
		mag(Rm, Gm); // Gm = |G|
		d1.fill_in(0);
		x1.fill_in(0);
		do 
		{
			x = x1;
			d = d1;
			x1 = Rm * x + z;
			mag(x1 - x, pom);
			d1 = Gm * d + pom + eps;
			mag(d1 - d, dpom);
		}
		while (!(dpom < eps));
		w = x0 + d1 * interval(-1.0, 1.0);
		return true;
	}
	return false;
}

bool Rohn(	const aafrmatrix& A, 
			const aafrvector& b, 
			ivector& w)
{
	int			n =	b.size();
	REAL		eps = 1.0e-16;
	dvector		x0(n), d(n), bm(n), d1(n), pom(n), dpom(n);
	dmatrix		R(n, n), M(n, n), I(n, n), Gm(n, n);
	ivector		x(n), x1(n), z(n);
    imatrix		C(n, n), Rm(n, n);
	AAFR		s, s2;
	//-------------------
	Rm = imatrix(n, n);
	mid(A, R);
	mid(b, bm);
	if (inv(R)) 
	{
		x0 = R * bm;
		for(int i = 0; i < n; i++) 
		{
			s2 = 0;
			for(int j = 0; j < n; j++) 
			{
				s = 0;
				for(int k = 0; k < n; k++)
					s = s + A(j, k) * x0[k];
				s2 = s2 + R(i, j) * (b[j] - s);
			}
			z[i] = s2.reduce();
		}
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++) 
			{
				s = 0;
				for(int k = 0; k < n; k++)
					s = s + R(i, k) * A(k, j);
				C(i, j) = s.reduce();
			}
		}
		Unit(I);
		Rm = I - C; // G = Rm
		mag(Rm, Gm);
		d1.fill_in(0);
		x1.fill_in(0);
		do 
		{
			x = x1;
			d = d1;
			x1 = Rm * x + z;
			mag(x1 - x, pom);
			d1 = Gm * d + pom + eps;
			mag(d1 - d, dpom);
		}
		while (!(dpom < eps));
		w = x0 + d1 * interval(-1.0, 1.0);
		return true;
	}
	return false;
}

bool RohnImpr(	const ivector& p, 
				const omatrix& oA, 
				const ovector& ob, 
				ivector& w2)
/* doesn't work properly and for now I don't know why */
{
	int			n =	ob.size(), np = p.size();
	REAL		eps = 1.0e-20;
	dvector		x0(n), d(n), d0(n), d1(n), d2(n), pom(n), pm(np), bm(n), dpom(n);
	dmatrix		R(n, n), I(n, n), Gm(n, n);
	ivector		x(n), x1(n), x2(n), z(n), w(n);
	imatrix		C(n, n), Rm(n, n);
	//-------------------
	w2 = ivector(n);
	mid(p, pm);
	bm = v_from_dep(pm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	if (inv(R)) 
	{
		x0 = R * bm; // to jest x*
		c_matrix(R, p, oA, C); // C = mid(A(p))^{-1}*A(p)
		z_vector(x0, R, p, oA, ob, z); // g = b = mid(A(p)^{-1}*b(q)
		Unit(I);
		Rm = I - C; // G = Rm
		mag(Rm, Gm);
		d1.fill_in(0);
		x1.fill_in(0);
		do 
		{
			x = x1;
			d = d1;
			x1 = Rm * x + z;
			mag(x1 - x, pom);
			d1 = Gm * d + pom + eps;
			mag(d1 - d, dpom);
		}
		while (!(dpom < eps));
		w = x0 + x + d * interval(-1.0, 1.0);
	}
	//=================
	x1 = w;
	d1 = d;
	x2 = Rm * x1 + z;
	d2 = Gm * d1;
	do 
	{
		x = x1; d = d1;
		x1 = x2; d1 = d2;
		x2 = Rm * x1 + z;
		d2 = Gm * d1;
		mag(d1 - d, dpom);
	}
	while (!(dpom < eps));
	w2 = x0 + x + d * interval(-1, 1);
	return true;
}

bool RohnImpr(	const ivector& p, 
				const ivector& q, 
				const omatrix& oA, 
				const ovector& ob, 
				ivector& w2)
/* doesn't work properly and for now I don't know why */
{
	int			n =	ob.size(), np = p.size(), nq = q.size();
	REAL		eps = 1.0e-20;
	dvector		x0(n), d(n), d0(n), d1(n), d2(n), pom(n), pm(np), qm(nq), bm(n), dpom(n);
	dmatrix		R(n, n), I(n, n), Gm(n, n);
	ivector		x(n), x1(n), x2(n), z(n), w(n);
	imatrix		C(n, n), Rm(n, n);
	//-------------------
	w2 = ivector(n);
	mid(p, pm);
	mid(q, qm);
	bm = v_from_dep(qm, ob); // midpoint vector
	R = m_from_dep(pm, oA); // midpoint matrix
	if (inv(R)) 
	{
		x0 = R * bm; // to jest x*
		c_matrix(R, p, oA, C); // C = mid(A(p))^{-1}*A(p)
		z_vector(x0, R, p, oA, ob, z); // g = b = mid(A(p)^{-1}*b(q)
		Unit(I);
		Rm = I - C; // G = Rm
		mag(Rm, Gm);
		d1.fill_in(0);
		x1.fill_in(0);
		do 
		{
			x = x1;
			d = d1;
			x1 = Rm * x + z;
			mag(x1 - x, pom);
			d1 = Gm * d + pom + eps;
			mag(d1 - d, dpom);
		}
		while (!(dpom < eps));
		w = x0 + x + d * interval(-1.0, 1.0);
	}
	x1 = w;
	d1 = d;
	x2 = Rm * x1 + z;
	d2 = Gm * d1;
	do 
	{
		x = x1; d = d1;
		x1 = x2; d1 = d2;
		x2 = Rm * x1 + z;
		d2 = Gm * d1;
		mag(d1 - d, dpom);
	}
	while (!(dpom < eps));
	w2 = x0 + x + d * interval(-1.0, 1.0);
	return true;
}