//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Interval base
//%
//%##########################################################################
//%

#include <iostream>
#include <cassert>
#include <cmath>
#include "interval_base.h"
#include "Round.h"

real nEta();
real nEps();
real vnEta = nEta();
real vnEps = nEps();

real
//%----------------------------------------------------------------------------
//% Calculates Epsilon (nEps = (min { x >= 0 : 1 + x > 1 }))
//%----------------------------------------------------------------------------
nEps()
{
	RoundNear();
	real c = 1, l, opc;
	do {
		l = c;
		c /= 2;
		opc = 1 + c;

	} while (opc > 1);
	return l;
}

real
//%----------------------------------------------------------------------------
//% Calculates Eta (nEta = min { x > 0 })
//%----------------------------------------------------------------------------
nEta()
{
	real c, l;
	c = nEps();
	do {
		l = c;
		c /= 2;
	} while (c > 0);
	return l;
}

real
//%----------------------------------------------------------------------------
//% Result: largest possible floating point number with r < a
//%----------------------------------------------------------------------------
predR(const real p)
{
	volatile real pr;
	RoundDown();
	pr = p - vnEta;
	RoundNear();
	return pr;
} //% >>> PredR <<<

real
//%----------------------------------------------------------------------------
//% Result: smallest possible floating point number with r > a
//%----------------------------------------------------------------------------
succR(const real p)
{
	volatile real pr;
	RoundUp();
	pr = p + vnEta;
	RoundNear();
	return pr;
} //% >>> SuccR <<<

static unsigned short
//%----------------------------------------------------------------------------
//% Definition of the smallest double
//%----------------------------------------------------------------------------
mask[16] = { /* bit masks for bits 0-15*/
	0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080,
	0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000
};

typedef union
{
	real dp; //% the 64-bit double precision value
	unsigned short sh[4]; //% overlay an array of short int
} xDouble;

#ifdef _MSC_VER
#define MSW 3	//% 3 if the left-most 16-bit short is sh[3]
#else
#define MSW 0   //% 0 if the left-most 16-bit short is sh[0]
#endif

real ulp(double x)
{
	xDouble	U, //% ulp for x
		X; //% working copy of X
	int	bit, //% position of bit e-1 in 16 word
		e1, //% biased exponent
		word; //% index of 16-bit word containing bit e-1
	X.dp = x;
	X.sh[MSW] &= 0x7FF0; //% isolate exponent in 16-bit word
	//% X.sh[0] now holds the exponent in bits 14-4
	U.dp = 0.0;
	if(X.sh[MSW] > 0x0340) { //% ulp is normalized number
		U.sh[MSW] = X.sh[MSW] - 0x0340; //% set exponent to e-52
	}
	//% the value 0x0340 is 52 left-shifted bits i.e. 0x0340=832=52<<4
	else { //% ulp is denormalized number
		e1 = (X.sh[MSW] >> 4) - 1; //% biased exponent - 1
		word = e1 >> 4; //% find 16-bit word containing bit e-1
		if(MSW == 0) word = 3 - word; //% compensate for word ordering
		bit = e1 % 16; //% find the bit position in the word
		U.sh[word] = mask[bit]; //% set the bit to 1
	}
	return U.dp; //% return ulp
}

union Double
{
	real x;
	struct {
#if ( ISLIB_ENDIANTEST == ISLIB_ONE )
		unsigned int mant1 : 32;
		unsigned int mant0 : 20;
		unsigned int expo : 11;
		unsigned int sign : 1;
#else
		unsigned int sign  : 1;
		unsigned int expo  : 11;
		unsigned int mant0 : 20;
		unsigned int mant1 : 32;
#endif
	} fx;
};

void testDoubleUnion()
{
	Double X;
	X.x = 2.0;
	std::cout << "sign = " << X.fx.sign << std::endl;
	std::cout << "exp = " << X.fx.expo << std::endl;
	std::cout << "mant0 = " << X.fx.mant0 << std::endl;
	std::cout << "mant1 = " << X.fx.mant1 << std::endl;
	X.x = sqrt(2.0);
	std::cout << "sign = " << X.fx.sign << std::endl;
	std::cout << "exp = " << X.fx.expo << std::endl;
	std::cout << "mant0 = " << X.fx.mant0 << std::endl;
	std::cout << "mant1 = " << X.fx.mant1 << std::endl;
	RoundUp();
	X.x = sqrt(2.0);
	std::cout << "sign = " << X.fx.sign << std::endl;
	std::cout << "exp = " << X.fx.expo << std::endl;
	std::cout << "mant0 = " << X.fx.mant0 << std::endl;
	std::cout << "mant1 = " << X.fx.mant1 << std::endl;
	RoundNear();
}

real
//%----------------------------------------------------------------------------
//% Smallest positive number
//%----------------------------------------------------------------------------
spn()
{
	Double X;
	X.fx.sign = 0;
	X.fx.expo = 1;
	X.fx.mant0 = 0;
	X.fx.mant1 = 0;
	return X.x;
}

real
//%----------------------------------------------------------------------------
//% 1ulp = 2^{e+1-N} , e = E - BIAS (BIAS = 1023), N = 53
//%----------------------------------------------------------------------------
ulpp(real x)
{
	const unsigned int BIAS = 1023;
	Double X;
	X.x = x;
	return ldexp(2.0, (X.fx.expo - BIAS) - 53);
}

real
//%----------------------------------------------------------------------------
//% 1ulp = 2^{e+1-N} , e = E - BIAS (BIAS = 1023), N = 53
//%----------------------------------------------------------------------------
ulp0(real x)
{
	Double X;
	unsigned int expx;

	X.x = x;
	expx = X.fx.expo;
	X.x = 0.0;
	std::cout << "exp = " << expx << std::endl;
	if (expx > 832) std::cout << "norm" << std::endl;
	X.fx.expo = (expx - 52);
	X.fx.mant1 = 0x0001;
	return X.x;
}

real
//%----------------------------------------------------------------------------
//% 1ulp = 2^{e-N} , e = E - BIAS (BIAS = 1023), N = 53
//%----------------------------------------------------------------------------
ulppp(real x)
{
	const unsigned int BIAS = 1023;
	Double X;
	X.x = x;
	return ldexp(2.0, (X.fx.expo - BIAS) - 52); // rather -52
}

int
//%----------------------------------------------------------------------------
//% http://en.wikipedia.org/wiki/Unit_in_the_last_place
//%----------------------------------------------------------------------------
nsig()
{
	float x = 1;
	int p = 0;
	while (x != x + 1) {
		x = x * 2;
		p = p + 1;
	}
	return p;
}

real
//%----------------------------------------------------------------------------
//% http://en.wikipedia.org/wiki/Unit_in_the_last_place
//%----------------------------------------------------------------------------
macheps()
{
	real machEps = 1.0;
	Double z;
	do {
		printf("%lf\t%.20lf\n", machEps, (1.0 + machEps));
		z.x = machEps;
		printf("%u %u %u %u\n", z.fx.sign, z.fx.expo, z.fx.mant0, z.fx.mant1);
		machEps /= 2.0;
		// If next epsilon yields 1, then break, because current
		// epsilon is the machine epsilon.
	} while ((1.0 + (machEps / 2.0)) != 1.0);
	return machEps;
}