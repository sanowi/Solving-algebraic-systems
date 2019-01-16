#pragma once
#ifndef _ROHN_H_
#define _ROHN_H_

class AAFR;

bool Rohn(const ivector&, const omatrix&, const ovector&, ivector&);
bool Rohn(const ivector&,	const ivector&, const omatrix&, const ovector&, ivector&);
bool Rohn(const aafrmatrix&, const aafrvector&, ivector&);
//
bool RohnImpr(const ivector&, const omatrix&, const ovector&, ivector&);
bool RohnImpr(const ivector&, const ivector&, const omatrix&, const ovector&, ivector&);

#endif