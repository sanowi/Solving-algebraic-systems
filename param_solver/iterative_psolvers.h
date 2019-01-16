#pragma once

class AAFR;

//%----------------------------------------------------------------------------
//% Auxiliary methods
//%----------------------------------------------------------------------------
int SpectralRadius(int op, int pr, const aafrmatrix& A0, dmatrix& B);
real SpectralRadius(int op, const imatrix& A);
void IterInner(const aafrvector& x0, dvector& innl, dvector& innu);
void IterInner(const aafrmatrix& x0, dmatrix& innl, dmatrix& innu);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems
//%----------------------------------------------------------------------------
int Rump(const aafrmatrix& A, const aafrvector& ob, ivector& w, ivector& innerw, int& niter);
int Rump(const aafrmatrix& A, const aafrmatrix& Bb, aafrmatrix& w, int& niter);
int Rump(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int RumpII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int Rump(const imatrix& A, const ivector& b, ivector& w);
int Rump(const ivector& p, const omvector& oA, const ovector& ob, ivector& w, int& stepCount);
//%----------------------------------------------------------------------------
int ParamGaussSeidel(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelRPrec(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelIII(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix R, aafrvector& w, int& niter);
int ParamGaussSeidel(const aafrmatrix& A, const aafrmatrix& B, aafrmatrix& x0, int& niter);
//int ParamGaussSeidelProba(const aafrmatrix& A, const aafrmatrix& B, aafrmatrix& x0, int& niter);
int MKolev(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevRC(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevRCII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
//%----------------------------------------------------------------------------
int Krawczyk(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int KrawczykII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int Krawczyk(const aafrmatrix& A, const aafrmatrix& bb, aafrmatrix& w, int& niter);
//%----------------------------------------------------------------------------
int SuccessiveOverrelaxation(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, real wgt, int& niter);