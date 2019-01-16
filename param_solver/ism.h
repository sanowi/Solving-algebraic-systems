#pragma once
#ifndef _ISM_H_
#define _ISM_H_

class AAFR;

bool Rohn(const ivector&, const omatrix&, const ovector&, ivector&);
bool Rohn(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool Rohn(const aafrmatrix&, const aafrvector&, ivector&);
bool RohnImpr(const ivector&, const omatrix&, const ovector&, ivector&);
bool RohnImpr(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);

bool HladikImpr(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool HladikImpr2(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool HladikImpr(const ivector&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool Hladik(const ivector&, const omatrix&, const ovector&, ivector&);
bool Hladik(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool Hladik(const aafrmatrix&, const aafrvector&, ivector&);

bool BauerSkeel(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool BauerSkeel(const aafrmatrix&, const aafrvector&, ivector&);
bool BauerSkeel(const ivector&, const omatrix&, const ovector&, ivector&);

bool HansenBliekRohn(const ivector&, const omatrix&, const ovector&, ivector&);
bool HansenBliekRohn(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool HansenBliekRohn(const aafrmatrix&, const aafrvector&, ivector&);

bool Proba(const ivector&, const omatrix&, const ovector&, ivector&);
bool Proba2(const ivector&, const omatrix&, const ovector&, ivector&);
bool Proba(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool Proba(const aafrmatrix&, const aafrvector&, ivector&);

bool ISM(const ivector&, const omatrix&, const ovector&, ivector&);
bool ISM(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&); /// ISM method with two sets of interval parameters
bool IISM(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&); /// 
bool ISM(const aafrmatrix&, const aafrvector&, ivector&);
bool ISM(const aafmatrix&, const aafvector&, ivector&);
bool ISM(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omatrix&, const ovector&, ivector&);
bool ISMModified(const ivector&, const omatrix&, const ovector&, ivector&);
bool ISM(const ivector&, const dmatrix&, const imatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool ISMC(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omatrix&, const ovector&, ivector&);

bool c_matrix(const dmatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const imatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const omatrix&, const ivector&, imatrix&);
bool b_matrix(const dmatrix&, const ivector&, const omatrix&, omatrix&);
bool b_vector(const dmatrix&, const ivector&, const ovector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool h_matrix(const dmatrix&);

#endif  // _ISM_H_
