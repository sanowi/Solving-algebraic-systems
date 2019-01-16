#include <iostream>
#include <cassert>
#include <cmath>
#include "config.h"
#include "gausselim.h"
#include "zvector.h"

bool c_matrix(const dmatrix& R, const ivector& p, const omatrix& oA, imatrix& C)
/* Calculates [C] = R * A([p]) */
{
	if (!oA.is_square()) return false;
    int			n = oA.num_rows();
    dvector		s(p.size());
	//-----------------------------
	for(int i = 0; i < n; i++) 
		for(int j = 0; j < n; j++) {
			s.fill_in((real)0);
			for(int k = 0; k < n; k++)
				s += oA(k, j) * R(i, k);
			if (s != (real)0) {
				C(i, j) = p * s;
			}
			else
				C(i, j) = (real)0;
		}
	return true;
}

bool c_matrix(const imatrix& R, const ivector& p, const omatrix& oA, imatrix& C)
/* Calculates [C] = [R] * A(p) */
{
	if (!oA.is_square()) return false;
    int			n = oA.num_rows();
    ivector		s(p.size());
	//-----------------------------
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++) {
			s.fill_in((real)0);
			for(int k = 0; k < n; k++)
				s += oA(k, j) * R(i, k);
			C(i, j) = p * s;
		}
	return true;
}

bool c_matrix(const omatrix& B, const ivector& p, imatrix& C)
/* Calculates [C] based on [B] = R * oA */
{
	if (!B.is_square()) return false;
	int n = B.num_rows();
	//-----------------------------
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			C(i, j) = p * B(i, j);
	return true;
}

bool b_matrix(const dmatrix& R, const ivector& p, const omatrix& oA, omatrix& B)
/* Calculates [B]  = R * oA */
{
	if (!oA.is_square()) return false;
    int			n = oA.num_rows(), np = p.size();
    dvector		s(np);
	//-----------------------------
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++) {
			s.fill_in((real)0);
			for(int k = 0; k < n; k++)
				s += oA(k, j) * R(i, k);
			B(i, j) = s;
		}
	return true;
}

bool z_vector(const dvector& x, const dmatrix& R, const ivector& p, 
			  const ivector& q, const omatrix& oA, const ovector& ob, ivector& z)
/*----------------------------------------*
 * Calculates [z] = R*([b(q)] - [A(p)]*x0)
 * R - real matrix, inverse of midpoint
 * p - left hand Vector of parameters
 * q - right hand Vector of parameters
 *----------------------------------------*/
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size(), nq = q.size();
    dvector		s1(nq), s2(np);
	//-----------------------------
	for(int i = 0; i < n; i++) {
		s1.fill_in((real)0);
		s2.fill_in((real)0);
		for(int j = 0; j < n; j++) {
			s1 += ob[ j ] * R(i, j);
			for(int k = 0; k < n; k++)
				s2 += oA(j, k) * (R(i, j) * x[ k ]);
		}
		z[ i ] = (q * s1) - (p * s2); // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
    }
	return true;
}

bool z_vector(const ivector& x, const imatrix& R, const ivector& p, 
			 const ivector& q, const omatrix& oA, const ovector& ob, ivector& z)
/*-------------------------------------------*
 * Calculates [z] = R*([b(q)] - [A(p)]*x0)
 * R - interval matrix, inverse of midpoint
 * p - left hand Vector of parameters
 * q - right hand Vector of parameters
 *-------------------------------------------*/
{
	if (!oA.is_square()) return false;
	int			n = ob.size(), np = p.size(), nq = q.size();
    ivector		s1(nq), s2(np);
	//-----------------------------
	for(int i = 0; i < n; i++) {
		s1.fill_in((real)0);
		s2.fill_in((real)0);
		for(int j = 0; j < n; j++) {
			s1 = s1 + ob[ j ] * R(i, j);
			for(int k = 0; k < n; k++)
				s2 = s2 + oA(j, k) * (R(i, j) * x[ k ]);
		}
		z[ i ] = q * s1 - p * s2; // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
    }
    return true;
}

bool z_vector(const dvector& b, const dvector& x, const dmatrix& R, 
			 const ivector& p, const omatrix& oA, const ovector& ob, ivector& z)
/*----------------------------------------*
 * Calculates [z] = R*([b(p)] - [A(p)]*x0)
 * p - left hand Vector of parameters
 *----------------------------------------*/
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size();
    dvector		s1(np), s2(np);
	//-----------------------------
	for(int i = 0; i < n; i++) {
         s1.fill_in((real)0);
		 s2.fill_in((real)0);
		 for(int j = 0; j < n; j++) {
              s1 += ob[ j ] * R(i, j);
              for(int k = 0; k < n; k++)
                   s2 += oA(j, k) * x[ k ] * R(i, j);
         }
		 z[ i ] = p * (s1 - s2);
    }
    return true;
}

bool z_vector(const dvector& x, const dmatrix& R, const ivector& p, 
			 const omatrix& oA, const ovector& ob, ivector& w)
{
	if (!oA.is_square()) return false;
	int			n = x.size(), np = p.size();
	real		s1, s2, s3, s4;
	interval	s5;
	//-----------------------------
	for(int i = 0; i < n; i++) {
		s1 = (real)0;
		s2 = (real)0;
		for(int j = 0; j < n; j++) {
			if(R(i, j) != (real)0 && ob[ j ][ 0 ] != (real)0)
				s1 += R(i, j) * ob[ j ][ 0 ];
			for(int k = 0; k < n; k++)
				s2 += R(i, j) * oA(j, k)[ 0 ] * x[ k ];
		}
		s1 -= s2;
		s5 = (real)0;
		for(int l = 1; l < np; l++) {
			s3 = (real)0;
			s4 = (real)0;
			for(int j = 0; j < n; j++) {
				s3 += R(i, j) * ob[ j ][ l ];
				for(int k = 0; k < n; k++)
					s4 += R(i, j) * oA(j, k)[ l ] * x[ k ];
			}
			s3 -= s4;
			s5 += s3 * p[ l ];
		}
		w[ i ] = s1 + s5;
	}
	return true;
}

bool h_matrix(const dmatrix& A)
/**-------------------------------------------------------------*
 * Verifies if a matrix A is an H-matrix
 * Def: A is an H-matrix if there exists a real vector u > 0, 
 * such that <A>u > 0, where <A> - Ostrovsky comparison matrix
 *--------------------------------------------------------------*/
{
	if (!A.is_square()) return false;
	int			n = A.num_rows();
	dvector		e(n), x(n);
	//-----------------------------
	e.fill_in((real)1);
	if (gw_gausse(A, e, x)) { // rozwiazanie ukladu <A>x = e
		for(int i = 0; i < n; i++) 
			if (x[ i ] <= 0) 
				return false; // not decided
		return true; // is H-matrix
	}
	return false; // not decided
}