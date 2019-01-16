/**--------------------------------------------------*
 *  Different versions of Skalna's Direct Method
 *---------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "config.h"
#include "aafr.h"
#include "oelement.h"
#include "inverse.h"
#include "gausselim.h"
#include "structure.h"
#include "ismmodf.h"

bool ISM(const ivector& p, const ivector& q, const omegamatrix& oA, const omegavector& ob, ivector& w)
/*-----------------------------------------------------------*
 * Skalna's Method for solving PILS
 * with two sets of interval parameters
 * - p - vector of left hand interval parameters
 * - q - vector of right hand interval parameters
 * - oA - matrix of left hand dependencies (A = oA * p)
 * - ob - vector of right hand dependencies (b = ob * p)
 *-----------------------------------------------------------*/
{
	if (!oA.is_square()) return false;
	int			n = oA.num_rows(), np = p.size(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq), pm(np);
	dmatrix		R(n, n), Cm(n, n);
    ivector		z(n);
    imatrix		C(n, n);
	//-----------------------------
	w = ivector(n);
	mid(q, qm);
	mid(p, pm);
	//time_t czas;
	//time(&czas);
	//std::cout << "przed vfromdep" << ctime(&czas) << std::endl;
	bm = v_from_dep(qm, ob);
	//time(&czas);
	//std::cout << "po vfromdep, przed mfromdep" << ctime(&czas) << std::endl;
	R = m_from_dep(pm, oA);
	//time(&czas);
	//std::cout << "po mfromdep, przed inv(R)" << ctime(&czas) << std::endl;
	if (inv(R)) {
		//time(&czas);
		//std::cout << "po inv(R), przed R*bm" << ctime(&czas) << std::endl;
		x = R * bm;
		//time(&czas);
		//std::cout << "po R*bm, przed cmatrix" << ctime(&czas) << std::endl;
		c_matrix(R, p, oA, C);
		//time(&czas);
		//std::cout << "po cmatrix, przed mig(C)" << ctime(&czas) << std::endl;
		mig(C, Cm);
		//time(&czas);
		//std::cout << "po mig(C), przed hmatrix(C)" << ctime(&czas) << std::endl;
		if (!h_matrix(Cm)) {/* H-matrix property verification */
			std::cout << "Probably not an H-matrix" << std::endl;
			return false; // not an H-matrix
		}
		//time(&czas);
		//std::cout << "po hmatrix(C), przed zvector" << ctime(&czas) << std::endl;
		z_vector(x, R, p, q, oA, ob, z); /* z = z(p, q) = R(b(q) - A(p)x0) */
		//time(&czas);
		//std::cout << "po zvector, przed mag(z)" << ctime(&czas) << std::endl;
		mag(z, magz);
		//time(&czas);
		//std::cout << "po mag(z), przed ge(Cm, magz)" << ctime(&czas) << std::endl;
		if (gw_gausse(Cm, magz, cmz)) {
			//time(&czas);
			//std::cout << "po ge(Cm, magz)" << ctime(&czas) << std::endl;
			for(int i = 0; i < n; i++)
				w[i] = x[i] + cmz[i] * interval(-1, 1); //interval(-cmz[i], cmz[i]);
			return true;
		}
	}
	std::cout << "Singular matrix!" << std::endl;
	return false; // singular matrix
}

//bool IISM(const ivector& p, const ivector& q, const omegamatrix& oA, const omegavector& ob, ivector& w)
///*--------------------------------------------------------------------*
// * Skalna's Method for solving parametric interval linear systems
// * with two Vectors of interval parameters
// * matrix R - the inverse of a midpoint matrix is calculated using 
// * the interval gaussian elimination
// *--------------------------------------------------------------------*/
//{
//	if (!oA.is_square()) return false;
//	int			n = oA.num_rows(), np = p.size(), nq = q.size();
//	dvector		bm(n), magz(n), qm(nq), pm(np);
//	dmatrix		Cm(n, n);
//	ivector		x(n), cmz(n), z(n), imagz(n);
//	imatrix		R(n, n), C(n, n), iCm(n, n);
//	//-----------------------------
//	w = ivector(n);
//	mid(q, qm);
//	mid(p, pm);
//	bm = v_from_dep(qm, ob);
//	to_int(m_from_dep(pm, oA), R);
//	if (inv(R)) {
//		x = R * bm;
//		c_matrix(R, p, oA, C);
//		mig(C, Cm);
//		if (!h_matrix(Cm)) /* H-matrix property verification */
//			return false; // not an H-matrix
//		z_vector(x, R, p, q, oA, ob, z);
//		mag(z, magz);
//		to_int(magz, imagz);
//		to_int(Cm, iCm);
//		if (gw_gausse(iCm, imagz, cmz)) {
//			for(int i = 0; i < n; i++)
//				w[i] = x[i] + cmz[i] * interval(-1, 1); //interval(-cmz[i], cmz[i]);
//			return true;
//		}
//	}
//	return false; // singular matrix
//}

bool ISM(const ivector& p, const ivector& q, const dmatrix& R, const imatrix& C, const omegamatrix& oA, const omegavector& ob, ivector& w)
/**----------------------------------------------------------*
 * ISM version with two sets of interval vectors
 * the difference is that R=(mid(A))^{-1} and C = <RA> 
 * are passed as arguments
 *-----------------------------------------------------------*/
{
	int			i, n = oA.num_rows(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq);
	dmatrix		Cm(n, n);
    ivector		z(n);
	//-----------------------------
	w = ivector(n);
	mid(q, qm);
	bm = v_from_dep(qm, ob);
	x = R * bm;
    mig(C, Cm);
	if (!h_matrix(Cm)) /* H-matrix property verification */
		return false; // not an H-matrix
    z_vector(x, R, p, q, oA, ob, z);
    mag(z, magz);
    if (gw_gausse(Cm, magz, cmz)) {
		for(i = 0; i < n; i++)
			w[i] = x[i] + cmz[i] * interval(-1, 1); //interval(-cmz[i], cmz[i]); //;
		return true;
    }
	return false; // singular matrix
}

bool ISMC(const ivector& p, const ivector& q, const dmatrix& R, const imatrix& C, const omegamatrix& oA, const omegavector& ob, ivector& w)
/*--------------------------------------------------------*
 * ISMF version with C = <RA> and R are being passed 
 * as an argument and with computing b = ob * q
 *--------------------------------------------------------*/
{
	if (!oA.is_square()) return false;
	int			i, n = oA.num_rows(), nq = q.size();
	dvector		x(n), bm(n), cmz(n), magz(n), qm(nq);
	dmatrix		Cm(n, n);
    ivector		b(n), z(n);
	//-----------------------------
	w = ivector(n);
	mid(q, qm);
	bm = v_from_dep(qm, ob); // tu zmiana 01.03.2009 - z pv na qv
	x = R * bm;
    mig(C, Cm);
	if (!h_matrix(Cm)) /* H-matrix property verification */
		return false; // not an H-matrix
    z_vector(x, R, p, q, oA, ob, z);
    mag(z, magz);
	if (gw_gausse(Cm, magz, cmz)) {
		for(i = 0; i < n; i++)
			w[i] = x[i] + cmz[i] * interval(-1, 1);
		return true;
	}
	return false; // singular matrix
}

//==============================================================
// zvector procedures
//==============================================================

bool c_matrix(const dmatrix& R, const ivector& p, const omegamatrix& oA, imatrix& C)
/* Calculates [C] = R * A([p]) */
{
	if (!oA.is_square()) return false;
    int		n = oA.num_rows();
    OEL		s;
	//-----------------------------
	for(int i = 0; i < n; i++) 
		for(int j = 0; j < n; j++) {
			s = 0;
			for(int k = 0; k < n; k++)
				s = s + oA(k, j) * R(i, k);
			C(i, j) = s * p;
		}
	return true;
}

//bool c_matrix(const imatrix& R, const ivector& p, const omegamatrix& oA, imatrix& C)
///* Calculates [C] = [R] * A(p) */
//{
//	if (!oA.is_square()) return false;
//    int		n = oA.num_rows();
//    OEL		s = 0;
//	//-----------------------------
//	for(int i = 0; i < n; i++)
//		for(int j = 0; j < n; j++) {
//			for(int k = 0; k < n; k++)
//				s = s + oA(k, j) * R(i, k);
//			C(i, j) = s * p;
//		}
//	return true;
//}

//bool c_matrix(const omegamatrix& B, const ivector& p, imatrix& C)
///* Calculates [C] based on [B] = R * oA */
//{
//	if (!B.is_square()) return false;
//	int n = B.num_rows();
//	//-----------------------------
//	for(int i = 0; i < n; i++)
//		for(int j = 0; j < n; j++)
//			C(i, j) = B(i, j) * p;
//	return true;
//}

bool b_matrix(const dmatrix& R, const ivector& p, const omegamatrix& oA, omegamatrix& B)
/* Calculates [B]  = R * oA */
{
	if (!oA.is_square()) return false;
    int		n = oA.num_rows(), np = p.size();
    OEL		s = 0;
	//-----------------------------
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++) {
			for(int k = 0; k < n; k++)
				s = s + oA(k, j) * R(i, k);
			B(i, j) = s;
		}
	return true;
}

bool z_vector(const dvector& x, const dmatrix& R, const ivector& p, 
			  const ivector& q, const omegamatrix& oA, const omegavector& ob, ivector& z)
/*----------------------------------------*
 * Calculates [z] = R*([b(q)] - [A(p)]*x0)
 * R - real matrix, inverse of midpoint
 * p - left hand Vector of parameters
 * q - right hand Vector of parameters
 *----------------------------------------*/
{
	if (!oA.is_square()) return false;
	int		n = ob.size(), np = p.size(), nq = q.size();
    OEL		s1, s2;
	//-----------------------------
	for(int i = 0; i < n; i++) {
		s1 = 0;
		s2 = 0;
		for(int j = 0; j < n; j++) {
			s1 = s1 + ob[j] * R(i, j);
			for(int k = 0; k < n; k++)
				if (R(i, j) != 0 && x[k] != 0)
					s2 = s2 + oA(j, k) * (R(i, j) * x[k]);
		}
		z[i] = s1 * q - s2 * p; // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
    }
	return true;
}

//bool z_vector(const ivector& x, const imatrix& R, const ivector& p, 
//			 const ivector& q, const omegamatrix& oA, const omegavector& ob, ivector& z)
///*-------------------------------------------*
// * Calculates [z] = R*([b(q)] - [A(p)]*x0)
// * R - interval matrix, inverse of midpoint
// * p - left hand Vector of parameters
// * q - right hand Vector of parameters
// *-------------------------------------------*/
//{
//	if (!oA.is_square()) return false;
//	int			n = ob.size(), np = p.size(), nq = q.size();
//    ivector		s1(nq), s2(np);
//	//-----------------------------
//	for(int i = 0; i < n; i++) {
//		s1.fill_in(0);
//		s2.fill_in(0);
//		for(int j = 0; j < n; j++) {
//			s1 = s1 + ob[j] * R(i, j);
//			for(int k = 0; k < n; k++)
//				s2 = s2 + oA(j, k) * (R(i, j) * x[k]);
//		}
//		z[i] = q * s1 - p * s2; // niezaleznie liczone jest dla wektora prawej strony i wektora lewej strony
//    }
//    return true;
//}

bool z_vector(const dvector& b, const dvector& x, const dmatrix& R, 
			 const ivector& p, const omegamatrix& oA, const omegavector& ob, ivector& z)
/*----------------------------------------*
 * Calculates [z] = R*([b(p)] - [A(p)]*x0)
 * p - left hand Vector of parameters
 *----------------------------------------*/
{
	if (!oA.is_square()) return false;
	int		n = oA.num_rows(), np = p.size();
    OEL		s1, s2;
	//-----------------------------
	for(int i = 0; i < n; i++) {
         s1 = 0;
		 s2 = 0;
		 for(int j = 0; j < n; j++) {
              s1 = s1 + ob[j] * R(i, j);
              for(int k = 0; k < n; k++)
                   s2 = s2 + oA(j, k) * x[k] * R(i, j);
         }
		 z[i] = (s1 - s2) * p;
    }
    return true;
}

//bool z_vector(const dvector& x, const dmatrix& R, const ivector& p, 
//			 const omegamatrix& oA, const omegavector& ob, ivector& w)
//{
//	if (!oA.is_square()) return false;
//	int			n = x.size(), np = p.size();
//	real		s1, s2, s3, s4;
//	interval	s5;
//	//-----------------------------
//	for(int i = 0; i < n; i++) {
//		s1 = 0;
//		s2 = 0;
//		for(int j = 0; j < n; j++) {
//			s1 = s1 + R(i, j) * ob[j][0];
//			for(int k = 0; k < n; k++)
//				s2 += R(i, j) * oA(j, k)[0] * x[k];
//		}
//		s1 -= s2;
//		s5 = 0;
//		for(int l = 1; l < np; l++) {
//			s3 = 0;
//			s4 = 0;
//			for(int j = 0; j < n; j++) {
//				s3 += R(i, j) * ob[j][l];
//				for(int k = 0; k < n; k++)
//					s4 += R(i, j) * oA(j, k)[l] * x[k];
//			}
//			s3 -= s4;
//			s5 += s3 * p[l];
//		}
//		w[i] = s1 + s5;
//	}
//	return true;
//}