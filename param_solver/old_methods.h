#pragma once

class AAFR;

bool NISM(const ivector&, const omatrix&, const ovector&, ivector&);
bool NISM_2(const ivector&, const omatrix&, const ovector&, ivector&);
bool NISM(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool NISM(const aafmatrix&, const aafvector&, ivector&, dmatrix& = dmatrix());
bool NISM(const aafrmatrix&, const aafrvector&, ivector&);

bool ISM(const ivector&, const omatrix&, const ovector&, ivector&);
bool ISM(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&); /// ISM method with two sets of interval parameters
bool IISM(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&); /// 

bool ISM(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omatrix&, const ovector&, ivector&);
bool ISMModified(const ivector&, const omatrix&, const ovector&, ivector&);
bool ISM(const ivector&, const dmatrix&, const imatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool ISMC(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omatrix&, const ovector&, ivector&);

bool ISM_I(const aafmatrix& A, const aafvector& b, ivector& w);
int ParamGaussSeidelNV(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int ParamGaussSeidelIINV(const aafrmatrix& A, const aafrvector& b, aafrvector& x0, int& niter);
int Rump(const aafrmatrix& A, const aafrvector& b, ivector& w, ivector& innerw, int& niter);
int Rump(const ivector& p, const omatrix& oA, const ovector& ob, ivector& w, int& stepCount);

bool c_matrix(const dmatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const dmatrix&, const ivector&, const omegamatrix&, imatrix&);
bool c_matrix(const imatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const omatrix&, const ivector&, imatrix&);
bool b_matrix(const dmatrix&, const ivector&, const omatrix&, omatrix&);
bool b_vector(const dmatrix&, const ivector&, const ovector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
int HMatrix(const dmatrix&);

bool ISMO(const ivector&, const omatrix&, const ovector&, ivector&);

