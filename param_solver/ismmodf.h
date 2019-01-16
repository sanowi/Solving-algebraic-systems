#pragma once
#ifndef _ISMMODF_H_
#define _ISMMODF_H_

class AAFR;

bool ISM(const ivector&, const ivector&, const omegamatrix&, const omegavector&, ivector&); /// ISM method with two sets of interval parameters
bool IISM(const ivector&, const ivector&, const omegamatrix&, const omegavector&, ivector&); /// 
bool ISM(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omegamatrix&, const omegavector&, ivector&);
bool ISM(const ivector&, const dmatrix&, const imatrix&, const ivector&, const omegamatrix&, const omegavector&, ivector&);
bool ISMC(const ivector&, const ivector&, const dmatrix&, const imatrix&, const omegamatrix&, const omegavector&, ivector&);

bool c_matrix(const dmatrix&, const ivector&, const omegamatrix&, imatrix&);
bool c_matrix(const imatrix&, const ivector&, const omegamatrix&, imatrix&);
bool c_matrix(const omatrix&, const ivector&, imatrix&);
bool b_matrix(const dmatrix&, const ivector&, const omegamatrix&, omegamatrix&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const ivector&, const omegamatrix&, const omegavector&, ivector&);
bool z_vector(const dvector&, const dvector&, const dmatrix&, const ivector&, const omegamatrix&, const omegavector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const omegamatrix&, const omegavector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const ivector&, const omegamatrix&, const omegavector&, ivector&);
bool h_matrix(const dmatrix&);

#endif  // _ISM_H_
