#pragma once

class AAFR;

//%----------------------------------------------------------------------------
//% Auxiliary methods
//%----------------------------------------------------------------------------
void IterInner(const aafrvector& x0, dvector& innl, dvector& innu);
void IterInner(const aafrmatrix& x0, dmatrix& innl, dmatrix& innu);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems
//%----------------------------------------------------------------------------
int Rump(const aafrmatrix& A, const aafrmatrix& Bb, aafrmatrix& w, int& niter);
int Rump(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int& niter);
int RumpRCond(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int& niter);
int RumpRLCond(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix& R, aafrvector& w, real &rho, int& niter);
int RumpRLCondLU(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int RumpRLCondLUP(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int RumpRLCondQR(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int RumpRLCondSVD(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int Rump(const aafrmatrix& A, const aafrmatrix& b, aafrmatrix& w, int& niter);
int Rump(const imatrix& A, const ivector& b, ivector& w);
int Rump(const ivector& p, const omvector& oA, const ovector& ob, ivector& w, int& stepCount);
//%----------------------------------------------------------------------------
int GaussSeidel(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidel(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelRC(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelRCond(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelRLCond(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix& R, aafrvector& w, int& niter);
int ParamGaussSeidel(const aafrmatrix& A, const aafrmatrix& B, aafrmatrix& x0, int& niter);
int PGSIRLCondLU(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int PGSIRLCondSVD(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
int PGSIRLCondQR(const aafrmatrix& A, const aafrvector& b, aafrvector& w, real &rho, int &niter);
//%----------------------------------------------------------------------------
int MKolev(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevRCond(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevRC(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int MKolevRCondRC(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
//%----------------------------------------------------------------------------
int Jacobi(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
//%----------------------------------------------------------------------------
int SOR(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, real wgt, int& niter);