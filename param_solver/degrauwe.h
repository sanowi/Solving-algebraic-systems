#pragma once

#include "../truss/TrussStructure.h"

int gem(int m, int op, const aafrmatrix& A, const aafrvector& b, aafrvector& w);
int gem(int m, int op, const aafrmatrix& A, const aafrvector& b, aafrvector& w, dmatrix& R);
int gem(int m, int op, const aafrmatrix& A, const aafrmatrix& b, aafrmatrix& w);
int gem(int m, const aafrmatrix& A, const aafrvector& b, aafrvector& w);
int gemMH(int m, const aafrmatrix& A, const aafrvector& b, aafrvector& w);
int gemMH(int mm, const aafrmatrix& A, const aafrmatrix& b, aafrmatrix& w);
void main_Degrauwe(Truss* ptruss, ivector& w, real& t);
void main_Degrauwe2Order(Truss* ptruss, ivector& w, real& t);
bool Degrauwe1Order(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, real& t);
bool Degrauwe2Order(const omvector& M1, const ovector& b1, const ivector& p1, ivector& w, real& t);