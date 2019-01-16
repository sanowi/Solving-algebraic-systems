#pragma once

class AAFR;

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems
//%----------------------------------------------------------------------------
int	BauerSkeel(const aafrmatrix& A, const aafrvector& b, ivector&  w);
int BauerSkeel(const dmatrix& R, const aafrmatrix& A, const aafrvector& b, ivector&  w);
int BauerSkeelRC(const aafrmatrix& A, const aafrvector& b, ivector&  w);
int BauerSkeelP(const aafrmatrix& A, const aafrvector& b, int r, ivector& w);
int BauerSkeel1(const aafrmatrix& A, const aafrvector& b, ivector&  w);
int BauerSkeel1RC(const aafrmatrix& A, const aafrvector& b, ivector& w);

//%----------------------------------------------------------------------------
//% Methods for solving parametric interval linear systems with multiple right-hand side
//%----------------------------------------------------------------------------
int BauerSkeel(const aafrmatrix& A, const aafrmatrix& B, imatrix& W);
int BauerSkeel(const aafrmatrix& A, const aafrmatrix& B, imatrix& W, dmatrix& R);

//%----------------------------------------------------------------------------
//% Methods for solving interval linear systems
//%----------------------------------------------------------------------------
int BauerSkeel(const imatrix& A, const ivector& b, ivector& w);
int BauerSkeelP(const imatrix& A, const ivector& b, ivector& w);
int BauerSkeelRC(const imatrix& A, const ivector& b, ivector& w);
int BauerSkeelPRC(const imatrix& A, const ivector& b, ivector& w);

