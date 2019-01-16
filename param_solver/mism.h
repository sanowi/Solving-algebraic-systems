#pragma once

class AAFR;

//%---------------------------------------------------
//% Auxiliary functions
//%---------------------------------------------------
bool h_matrix(const dmatrix&);
void OUSystem(const aafrmatrix& A, const aafrvector& b, aafrmatrix& AA, aafrvector& bb);

//%---------------------------------------------------
//% Methods
//%---------------------------------------------------
bool ISM(const aafrmatrix&, const aafrvector&, ivector&);
bool ISMP0(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
bool ISMP(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
bool ISM(const aafrmatrix& A, const aafrvector& b, ivector& w, dmatrix& R);
bool ISM(const aafrmatrix& A, const aafrmatrix& B, imatrix& w);
bool ISM(const aafrmatrix& A, const aafrmatrix& B, imatrix& w, dmatrix& R);



