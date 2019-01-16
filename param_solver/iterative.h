#pragma once

class AAFR;

bool spectralRadius(int op, const aafrmatrix& A0, dmatrix& B);
void iterInner(const aafrvector& x0, dvector& innl, dvector& innu);
void iterInner(const aafrmatrix& x0, dmatrix& innl, dmatrix& innu);

bool gsOuter(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
bool gsOuterII(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
bool gsOuter(const aafrmatrix& A, const aafrmatrix& B, aafrmatrix& x0, int& niter);

bool mkolev(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool mkolevII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool mkolevRC(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool mkolevRCII(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);

bool Krawczyk(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool Krawczyk2(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool Krawczyk(const aafrmatrix& A, const aafrmatrix& bb, aafrmatrix& w, int& niter);

bool jacobi(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
bool jacobi(const imatrix& A, const ivector& b, ivector& w);
bool sor(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, real wgt, int& niter);
bool RohnOverdeterm(const aafrmatrix& A, const aafrvector& b, ivector& w);