#pragma once

class AAFR;

//%----------------------------------------------------------------------------
//% Methods for solving interval linear systems
//%----------------------------------------------------------------------------
int FirstSolution(const dmatrix& R, const imatrix& A, const ivector& b, ivector& w);
//%----------------------------------------------------------------------------
int GaussSeidel(const imatrix&, const ivector&, ivector&, int& niter);
int GaussSeidelRC(const imatrix&, const ivector&, ivector&, int& niter);
//%----------------------------------------------------------------------------
int Jacobi(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int Jacobi(const imatrix& A, const ivector& b, ivector& w, int& niter);
int JacobiRC(const imatrix& A, const ivector& b, ivector& w, int& niter);
//%----------------------------------------------------------------------------
int Krawczyk(const imatrix&A, const ivector& b, ivector& w, int& niter);
int KrawczykRC(const imatrix&A, const ivector& b, ivector& w, int& niter);