#ifndef _SensAnal_H_
#define _SensAnal_H_

ivector SAnal(const ivector&, const omatrix&, const ovector&);
ivector SAnalF(const ivector&, const ivector&, const omatrix&, const ovector&);
void SAnalMonotTest(const ivector& pv, const ivector& qv, const omatrix& oA, const ovector& ob);
ivector SAnalFW(const ivector&, const ivector&, const omatrix&, const ovector&);

#endif