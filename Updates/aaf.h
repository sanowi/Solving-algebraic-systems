//%
//%=========================================================
//%
//% Affine Forms (AAF) Class
//% Written by: Iwona Skalna
//%
//% An affine form is a degree-1 polynomial:
//% x = x_0 + x_1*e_1 + ... + x_n*e_n
//% where:
//% n - the length of an affine form
//% x_0 - central value
//% x_i - partial deviations, i = 1, ..., n
//% e_i - noise symbols (e_i in [-1, 1]), i = 1, ..., n
//%
//% Remark: affine form is defined by two vectors
//% vector of indices (in increasing order) and vector of coefficients
//% nly non-zero elements are stored
//% (the idea resembles sparse vectors)
//%
//%----------------------------------------------------------
//%

#pragma once

enum MULTS_CASES {
	MULTS_QUICK = 1,//% QUICK MULT
	MULTS_BEST = 2, //% BEST MULT
	MULTS_STD = 3,	//% STD MULT
	MULT_MINR = 4,	//% MIN RANGE
	MULTS_MH = 5,	//% Hladik's approach
	MULTS_MHS = 6,	//% Hladik's approach with scaling vector
	MULTS_SAA = 7,	//% Trade of between accuracy and complexity
	MULTS_AASEE = 8 //% Tangent plane Best??? (for some problems maybe)
};

class AAF
{
protected:
	dvector	_coeffs;			//% partial deviations
	cvector	_indexes;			//% only non-zero indices are stored
	unsigned int _length;		//% length of an affine form
	static unsigned int	_last;	//% number of the highest noise symbol in use
public:
	AAF(real c = 0.0) //% AAF::AAF( double c = 0.0 )
	{
		//% constructor 1; creates an affine form from a real number
		_coeffs = dvector(1);
		_indexes = cvector(1);
		_coeffs[0] = c;
		_indexes[0] = 0;
	}
	AAF(cvector idx, dvector cfs) : _indexes(idx), _coeffs(cfs)
	{
		//% constructor 2; creates an affine form from a vector 
		//% of indices and a vector of coefficients
	}
	AAF(const interval& a); //% constructor 3; creates an affine form from an interval
	AAF(const AAF&); //% constructor 4; copying constructor 
	virtual ~AAF() {
		//%  destructor
	}
	//%------------------------------------------------------------------------------------
	static cvector idx_union(cvector&, cvector&);		//% set union of indices
	static dvector cfs_union(AAF&, AAF&);				//% formerly set union of coefficients
	static int inclast() { return ++_last; };			//% increments _last used index
	static void setlast(int val) { _last = val; };		//% sets _last use index
	static int getlast() { return _last; };				//% gets _last used index
	static AAF quickMult(const AAF&, const AAF&);		//% trivial multiplication of affine forms
	static AAF bestMult(const AAF&, const AAF&);		//% Chebyshev multiplication of affine forms
	static AAF stdMult(const AAF&, const AAF&);			//% standard (Kolev) multiplication of affine forms
	static AAF stdMultImp(const AAF&, const AAF&);		//% standard (Kolev) multiplication of affine forms
	static AAF hlMult(const AAF&, const AAF&);			//% Hladik's multiplication
	static AAF hlMultScaled(const AAF&, const AAF&);	//% Hladik's multiplication with scaling vector
	static AAF AASEEMult(const AAF&, const AAF&);
	//%------------------------------------------------------------------------------------
	AAF& operator = (const AAF&);						//% copying operator
	friend AAF operator + (const AAF&, const AAF& b);
	friend AAF operator + (const AAF&, const real);
	friend AAF operator + (const real, const AAF&);
	friend AAF operator - (const AAF&);
	friend AAF operator - (const AAF&, const AAF&);
	friend AAF operator - (const AAF&, const real);
	friend AAF operator - (const real, const AAF&);
	friend AAF operator * (const AAF& a, const AAF&);
	friend AAF operator * (const real, const AAF&);
	friend AAF operator * (const AAF&, const real);
	friend AAF operator / (const AAF&, const AAF&);
	friend void	operator += (AAF&, const AAF&);
	friend void	operator += (AAF&, const real);
	friend std::ostream& operator << (std::ostream& os, const AAF&);
	//%------------------------------------------------------------------------------------
	cvector	idx() const { return _indexes; } //% returns indices
	dvector	cfs() const { return _coeffs; }  //% returns coefficients
	real mid() { return _coeffs[0]; }		 //% midpoint of an affine form
	real rad();				 //% radius of an affine form
	interval reduce() const; //% converts affine form to interval
	void reduceZeroes();	 //% reduces coefficients equal to zero
	void setZeroAt(int);	 //% set i-th coefficient to zero

	AAF	reciprocal() const;  //% reciprocal of affine form - Chebyshev
	AAF	reciprocal2() const; //% reciprocal of affine form - min-range
	AAF	sqrtf();			 //% square root of affine form - Chebyshev
	AAF	sqrttf();			 //% square root of affine form - min-range
	AAF	sqrf();				 //% square of affine form - Chebyshev
	AAF	sqrrf();			 //% square of affine form - min-range
	AAF	expf();				 //% exponent of affine form - Chebyshev
}; //% AAF class
