#pragma once
#ifndef _IMLLIB_H_
#define _IMLLIB_H_

int CG(const dmatrix &A, dvector &x, const dvector &b, dmatrix &M, int &max_iter, real &tol);
int CG(const aafmatrix &A, ivector &x, const aafvector &b, dmatrix &M, int &max_iter, real &tol);

#endif
