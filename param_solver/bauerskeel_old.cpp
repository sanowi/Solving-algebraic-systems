#include "stdafx.h"
#include "config.h"
#include "aafr.h"
#include "inverse.h"
#include "mism.h"
#include "gausselim.h"
#include "TrussStructure.h"
#include "bauerskeel.h"

// czesciowa specjalizacja szanlonu (template)

bool BauerSkeel(
			const ivector& p, 
			const omatrix& oA, 
			const ovector& ob, 
			ivector& w)
{
	int			n =	ob.size(), np = p.size();
	dvector		pm(np), bm(n), br(n), xc(n), y(n), w0(n), wup(n), wdown(n), mxc(n);
	dmatrix		R(n, n), M(n, n), Rm(n, n), I(n, n);
	ivector		b(n);
    imatrix		C(n, n);
	//-------------------
	mid(p, pm);
	bm = v_from_dep(pm, ob);		// bm = b(mid(p))
	R = m_from_dep(pm, oA);			// R = A(mid(p))
	if (inv(R))						// R = A(mid(p))^{-1}
	{
		xc = R * bm;				// midpoint solution: x_c= mid(A(p))^{-1}*mib(b)
		mag(xc, mxc);				// mxc = |x_c|
		// preconditioning
		c_matrix(R, p, oA, C);		// C = R*A(p)
		b_vector(R, p, ob, b);		// b = R*b(p)
		rad(b, br);					// br = rad(b)
		rad(C, M);					// M = rad(C) = |R|*rad(A)
		Unit(I);
		Rm = I - M;					// Rm = I-|R|*rad(A)
		if (inv(Rm))				// now Rm = (I - |R|*rad(A))^{-1}
		{ 
			y = Rm * (mxc + br);
			wdown = xc - y + mxc;
			wup = xc + y - mxc;
			for(int i = 0; i < n; i++)
			{
				w[i] = interval(wdown[i], wup[i]);
			}
			return true;
		}
	}
	return false;
}

bool BauerSkeel(const aafrmatrix& A, const aafrvector& B, ivector& w)
{
	int			n =	B.size();
	dvector		bm(n), br(n), x0(n), y(n), w0(n), wup(n), wdown(n), mx0(n);
	dmatrix		R(n, n), M(n, n), Rm(n, n), I(n, n);
	ivector		b(n);
    imatrix		C(n, n);
	AAFR		s;
	//-------------------
	w = ivector(n);
	mid(A, R);
	mid(B, bm);
	if (inv(R)) 
	{
		x0 = R * bm;		// x0 = mid(A(p))^{-1}*mib(b)
		mag(x0, mx0);		// to jest |x*|
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++) 
			{
				s = 0.0;
				for(int k = 0; k < n; k++)
				{
					s = s + R(i, k) * A(k, j);
				}
				C(i, j) = s.reduce();
			}
		}
		for(int i = 0; i < n; i++) 
		{
			s = 0;
			for(int j = 0; j < n; j++)
			{
				s = s + R(i, j) * B[j];
			}
			b[i] = s.reduce();
		}
		rad(b, br); // br = rad(b)
		rad(C, M); // M = rad(C)
		Unit(I);
		Rm = I - M;
		if (inv(Rm)) // Rm = (I - M)^{-1}
		{ 
			y = Rm * (mx0 + br);
			wdown = x0 - y + mx0;
			wup = x0 + y - mx0;
			for(int i = 0; i < n; i++)
			{
				w[i] = interval(wdown[i], wup[i]);
			}
			return true;
		}
	}
	return false;
}

ivector BauerSkeel(	const aafrmatrix& A, const aafrvector& b)
/*------------------------------------------------*
 * BauerSkeel version with computations performed
 * directly on affine forms
 *------------------------------------------------*/
{
	if (!A.is_square()) return false;
	int n = A.num_rows();
	ivector		w(n), z(n);
	imatrix		C(n, n);
	dmatrix		R(n, n), Cm(n, n);
	dvector		bm(n), x(n), magz(n), cmz(n);
	AAFR		s, s2;
	//-----------------------------
	w = ivector(n);
	w.fill_in(0.0);
	mid(A, R);
	mid(b, bm);
	if (inv(R)) 
	{
		x = R * bm;
		for(int i = 0; i < n; i++) 
		{
			s2 = 0.0;
			for(int j = 0; j < n; j++) 
			{
				s = 0.0;
				for(int k = 0; k < n; k++)
				{
					s = s + A(j, k) * x[k];
				}
				s2 = s2 + R(i, j) * (b[j] - s);
			}
			z[i] = s2.reduce();
		}
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++) 
			{
				s = 0.0;
				for(int k = 0; k < n; k++)
				{
					s = s + R(i, k) * A(k, j);
				}
				C(i, j) = s.reduce();
			}
		}
		mig(C, Cm);
		if (!h_matrix(Cm)) /* H-matrix property verification */
		{
			return false; // not an H-matrix
		}
		mag(z, magz);
		if (gw_gausse(Cm, magz, cmz)) 
		{
			for(int i = 0; i < n; i++)
			{
				w[i] = x[i] + cmz[i] * interval(-1, 1);
			}
			return true;
		}
	}
	return false; // singular matrix
}