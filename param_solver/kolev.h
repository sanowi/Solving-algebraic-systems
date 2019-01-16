#pragma once

void reduce(const dvector& t1, const dmatrix& t2, const dvector& t3, ivector& v);

void eIntPLP(int minmax, const dvector& c, const dvector& x0, const dvector& t1, const dmatrix& t2, const dvector& t3, interval& eI);

void kolevReducei(int i, const dvector& t1, const dmatrix& t2, const dvector& t3, interval& v);

void omvec2real(const omvector& M, const dvector& p, dmatrix& A);
void ovec2real(const ovector& b, const dvector& p, dvector& bb);

//void kolevorg(const omvector& M, const ovector& b, const ivector& p, ivector& w, dvector& t1, dvector& t3, dmatrix& t2);
//void kolev(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dvector& x0);

void kolevSWIM(int minmax, const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, int fun);
bool kolevOuter(const omvector& M, const ovector& b, const ivector& p, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dvector& x0, int& niter);
bool kolevOuter(int mode, const omvector& M1, const ovector& b1, const ivector& p1, ivector& w);
bool kolevOuter(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w);
bool kolevOuter(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dvector& x0, int& niter);

void kolevInner(const dvector& x0, const dvector& t1, const dvector& t3, const dmatrix& t2, ivector& innsol);
void kolevInner(const dvector& x0, const dvector& t1, const dvector& t3, const dmatrix& t2, dvector& innl, dvector& innu);

bool kolevDirect(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dvector& x0);

void kolevL(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dvector& x0); // Sergio
void kolevQ(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w);
void kolevQ(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, dvector& t1, dvector& t3, dmatrix& t2, dmatrix& t4, dvector& x0);

void ktransform(const omvector& M1, const ovector& b1, const ivector& p1, omvector& M, ovector& b, ivector& p);
void ktransform(const omvector& M1, const ivector& p1, omvector& M, ivector& p);
void ktransformEig(const omvector& M1, const ivector& p1, omvector& M, ivector& p);
real stopCrit(const ivector& x, const ivector& y);
