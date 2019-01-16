#pragma once

class AAFR;

bool BauerSkeel(const aafrmatrix& A, const aafrvector& b, int r, ivector&  w);
bool BauerSkeelP(const aafrmatrix& A, const aafrvector& b, int r, ivector& w);
bool BauerSkeel1(const aafrmatrix& A, const aafrvector& b, ivector&  w);
bool BauerSkeel1RC(const aafrmatrix& A, const aafrvector& b, ivector& w);

bool BauerSkeel(const aafrmatrix& A, const aafrmatrix& B, imatrix& W);
bool BauerSkeel(const aafrmatrix& A, const aafrmatrix& B, imatrix& W, dmatrix& R);

bool BauerSkeel(const ivector&, const omatrix&, const ovector&, ivector&);
bool BauerSkeel(const imatrix& A, const ivector& b, ivector& w);
bool BauerSkeelP(const imatrix& A, const ivector& b, ivector& w);
bool BauerSkeelRC(const imatrix& A, const ivector& b, ivector& w);

bool Hladik(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool Hladik(const aafmatrix&, const aafvector&, ivector&);
bool HladikImpr(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool HladikImpr2(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool HladikImpr(const ivector&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
