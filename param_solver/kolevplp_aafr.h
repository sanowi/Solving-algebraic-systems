#pragma once

void kolevPLPOuter(const dvector& c, const aafrmatrix& M1, const aafrvector& b1, interval& mw, real& time);

void kolevPLP(int minmax, const dvector& c, const aafrmatrix& M1, const aafrvector& b1, interval& mw, int& iter, real& time);