//%
//%==========================================================
//%
//% Revised Affine Forms (RAF) class
//% an affine form is a first-order polynomial
//% x = x_0 + x_1*e_1 + ... + x_n*e_n + r
//% - n - constant length of the form
//% - x_0 - central value
//% - x_i - partial deviations
//% - e_i - dummy variables (e_i in [-1, 1])
//%
//%==========================================================
//%

#pragma once

#include "../vector_matrix/vector_matrix_base.h"
#include "aaf.h"

enum MULT_CASES {
	MULT_BEST = 1,
	MULT_NEWBEST = 2,
	MULT_STD = 3,
	MULT_QUICK = 4,
	MULT_MIY = 5,
	MULT_RUKAS = 6,
	MULT_NEWBESTON = 7,
	MULT_AASEE = 8
};

//% 1 - Chebyshev approximation (CA)
//% 2 - New Chebyshev multiplication (NCA)
//% 3 - Standard multiplication (SM)
//% 4 - Trivial multiplication (TM)
//% 5 - Modification of SM by Miyajima
//% 6 - Modification Rump-Kashiwagi
//% 7 - spectral decomposition

class AAFR : public AAF
{
private:
	real _r; //% accumulative error
public:
	AAFR(real c = 0) : AAF(c) { _r = 0.0; } //% constructor 1, from a real number
	AAFR(cvector _indexes, dvector _coeffs, real r = 0.0) : 
		AAF(_indexes, _coeffs) { _r = r; } //% constructor 2, from indexes and coefficients
	AAFR(interval& a) : AAF(a) { _r = 0.0; } //% constructor 3, from an interval
	AAFR(real c, real r) : AAF(c) { _r = r; } //% constructor 4, from two real numbers
	AAFR(const AAFR&);
	virtual ~AAFR(){} //% destructor
	AAFR& operator = (const AAFR&);
	interval reduce() const;
	AAFR sqrtf();
	AAFR reciprocal() const;
	AAFR reciprocal_mr() const;
	AAFR reciprocal2();
	AAFR sqr();
	AAFR cube();
	AAFR quickmult(const AAFR&);
	real r() const { return _r; }
	void blow(real mi, real e);
	//%--------------------------------------------------------
	static AAFR bestMult(const AAFR&, const AAFR&);		 //% Chebyshev multiplication of revised affine forms
	static AAFR newBestMult(const AAFR&, const AAFR&);	 //% O(n*log(n)) slope version of Chebyshev multiplication of revised affine forms
	static AAFR newBestMultON(const AAFR&, const AAFR&); //% O(n) slope version of Chebyshev multiplication of revised affine forms
	static AAFR stdMult(const AAFR&, const AAFR&);		 //% standard (Kolev) multiplication of affine forms
	static AAFR quickMult(const AAFR&, const AAFR&);	 //% trivial multiplication of affine forms
	static AAFR miyMult(const AAFR&, const AAFR&);		 //% Hladik's multiplication
	static AAFR rukasMult(const AAFR&, const AAFR&);	 //% Rump-Kashiwagi multiplication
	static AAFR quickMultKolev(const AAFR&, const AAFR&);//% Standard multiplication
	static AAFR AASEMult(const AAFR&, const AAFR&);
	static AAFR AASEMult2(const AAFR&, const AAFR&);
	//%--------------------------------------------------------
	friend AAFR operator + (const AAFR& a, const AAFR& b);
	friend AAFR operator + (const AAFR& a, const real b);
	friend AAFR operator + (const real a, const AAFR& b);
	friend AAFR operator + (const AAFR& a, const interval b);
	friend AAFR operator - (const AAFR& a);
	friend AAFR operator - (const AAFR& a, const AAFR& b);
	friend AAFR operator - (const AAFR& a, const real b);
	friend AAFR operator - (const real a, const AAFR& b);
	friend AAFR operator * (const AAFR& a, const AAFR& b);
	friend AAFR operator * (const real a, const AAFR& b);
	friend AAFR operator * (const AAFR& a, const real b);
	friend AAFR operator / (const AAFR& a, const AAFR& b);
	friend void operator += (AAFR& a, const AAFR& b);
	friend void operator += (AAFR& a, const real b);
	friend bool operator == (const AAFR& a, const real b);
	friend std::ostream& operator << (std::ostream& os, const AAFR&);
}; //% >>> AAFR class <<<
