#pragma once

class AAFR;

//%////////////////////////////////////////////////////////////////////////////
//% Auxiliary functions
//%----------------------------------------------------------------------------
int HMatrix(const dmatrix&);
void OUSystem(const aafrmatrix& A, const aafrvector& b, aafrmatrix& AA, aafrvector& bb);

//%////////////////////////////////////////////////////////////////////////////
//% Methods
//%----------------------------------------------------------------------------
int ParametricDirectMethod(const aafrmatrix&, const aafrvector&, ivector&);
int ParametricDirectMethodPII(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
int ParametricDirectMethodPIII(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix& R, ivector& w);
int ParametricDirectMethodPIII2(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix& R, ivector& w);
int ParametricDirectMethodPP(const aafrmatrix& A, const aafrvector& b, const dmatrix& L, const dmatrix& R, ivector& w);
int ParametricDirectMethodP(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
int ParametricDirectMethod(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
int ParametricDirectMethod(const aafrmatrix& A, const aafrmatrix& B, imatrix& w);
int ParametricDirectMethod(const aafrmatrix& A, const aafrmatrix& B, imatrix& w, dmatrix& R);
int ParametricDirectMethod(const imatrix& A, const ivector& b, ivector& w);

int KParametricDirectMethod(const aafrmatrix& A, const aafrvector& b, aafrvector& w);
int KParametricDirectMethodII(const aafrmatrix& A, const aafrvector& b, aafrvector& w);
int KParametricDirectMethod(const aafrmatrix& A, const aafrmatrix& b, aafrmatrix& w);



