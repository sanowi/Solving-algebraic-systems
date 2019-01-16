#pragma once

class AAFR;

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems
//%----------------------------------------------------------------------------
int HansenBliekRohn1(const aafrmatrix& A, const aafrvector& b, ivector& result);
int HansenBliekRohn1RC(const aafrmatrix& A, const aafrvector& b, ivector& result);
int HansenBliekRohn(const aafrmatrix& A, const aafrvector& b, ivector& result); 
int HansenBliekRohn(const aafrmatrix& A, const aafrvector& b, ivector& result, dmatrix& midinv);
int HansenBliekRohnKRI(const aafrmatrix& A, const aafrvector& b, ivector& result, aafrmatrix& V, aafrvector& v, dmatrix& midinv);
int HansenBliekRohnRC(const aafrmatrix& A, const aafrvector& b, ivector& result);
int HansenBliekRohnP(const aafrmatrix& A, const aafrvector& b, int r, ivector& result);
int HansenBliekRohnComb(const aafrmatrix& A, const aafrvector& b, ivector& result);
//%----------------------------------------------------------------------------
int NingKearfott(const aafrmatrix& A, const aafrvector& b, ivector& w);
int NingKearfottRC(const aafrmatrix& A, const aafrvector& b, ivector& w);
int NingKearfott(const imatrix& A, const ivector& b, ivector& w);
//%----------------------------------------------------------------------------
int NewOperatorMH(const aafrmatrix& A, const aafrvector& b, ivector& w);
int NewOperatorMHRC(const aafrmatrix& A, const aafrvector& b, ivector& w);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems
//% that used right preconditioning
//%----------------------------------------------------------------------------
int HansenBliekRohnII0(const aafrmatrix& A, const aafrvector& b, ivector& result, dmatrix& midinv);
int HansenBliekRohnII(const aafrmatrix& A, const aafrvector& b, ivector& result);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems with multiple right-hand side
//%----------------------------------------------------------------------------
int HansenBliekRohnMRHS(const aafrmatrix& A, const aafrmatrix& b, imatrix& result);

//%----------------------------------------------------------------------------
//% Auxiliary methods
//%----------------------------------------------------------------------------
void BS_vs_HBR(const dmatrix& M, const dvector& xc, const dvector& xa, dvector& errl, dvector& erru);

//%----------------------------------------------------------------------------
//% Methods for solving interval linear systems
//%----------------------------------------------------------------------------
int HansenBliekRohn(const imatrix& A, const ivector& b, ivector& result, int op);
int HansenBliekRohnRC(const imatrix& A, const ivector& b, ivector& result);
int HansenBliekRohnP(const imatrix& A, const ivector& b, ivector& result);
int HansenBliekRohnRCP(const imatrix& A, const ivector& b, ivector& result);
//%----------------------------------------------------------------------------
int NingKearfottRC(const imatrix& A, const ivector& b, ivector& w);
//%----------------------------------------------------------------------------
int NewOperatorMH(const imatrix& A, const ivector& b, ivector& w);
int NewOperatorMHRC(const imatrix& A, const ivector& b, ivector& w);
//%----------------------------------------------------------------------------
int ErrorImpr(const aafrmatrix& A, const aafrvector& b, ivector& w);
