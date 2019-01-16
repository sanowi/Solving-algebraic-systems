#pragma once
#ifndef _RUMP_FPI_H_
#define _RUMP_FPI_H_

bool RumpInner(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool RumpOuter(const ivector&, const omatrix&, const ovector&, ivector&);
bool RumpOuter(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool RumpOuterGSS(const ivector&, const omatrix&, const ovector&, ivector&);
bool RumpOuterGSS(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
bool RUMPInnerF(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);
//ivector RumpOuter(const imatrix&, const ivector&);
//dvector RumpPopova(const ivector&, const omatrix&, const ovector&, int, int);

#endif