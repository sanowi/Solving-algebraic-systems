#pragma once

class AAFR;

int Rump(const aafrmatrix& A, const aafrvector& ob, ivector& w, ivector& innerw, int& niter);
int Rump(const aafrmatrix& A, const aafrmatrix& Bb, aafrmatrix& w, int& niter);
int Rump(const aafrmatrix& A, const aafrvector& b, aafrvector& w, int& niter);
int Rump(const imatrix& A, const ivector& b, ivector& w);
