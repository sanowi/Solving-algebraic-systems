#pragma once
#ifndef _ZVECTOR_H_
#define _ZVECTOR_H_

bool c_matrix(const dmatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const imatrix&, const ivector&, const omatrix&, imatrix&);
bool c_matrix(const omatrix&, const ivector&, imatrix&);
bool b_matrix(const dmatrix&, const ivector&, const omatrix&, omatrix&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const dvector&, const dmatrix&, const ivector&, const omatrix&, const ovector&, ivector&);
bool z_vector(const ivector&, const imatrix&, const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool h_matrix(const dmatrix&);

#endif  // _zVector_H_
