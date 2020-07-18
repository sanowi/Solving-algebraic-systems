//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Module: Elementary interval functions
//%
//%##########################################################################
//%

#include <iostream>
#include <fstream>
#include <conio.h>
#include <math.h>
#include "../interval/interval_base.h"
#include "../interval/interval.h"
#include "../vector_matrix/vector_matrix_base.h"
#include "intfun.h"

using namespace islib_interval;

static int
//%----------------------------------------------------------------------------
//% Calculates the quadrant x is lying in (x in [0, 2Pi])
//%----------------------------------------------------------------------------
Quadrant(const real x)
{
	if (x <= hPi) return 0; // half pi
	if (x <= Pi) return 1; // pi
	if (x <= Pi + hPi) return 2; // one and half pi
	return 3; // greater than one and half pi
}

real
//%----------------------------------------------------------------------------
//% Moves the angle 'a' but in the interval [0,2pi]
//%----------------------------------------------------------------------------
NormalizeAngle(real new_a)
{
	//REAL new_a = a;

	if ((new_a > Pi) || (new_a < -Pi))
		new_a = atan2(sin(new_a), cos(new_a));

	if (new_a < 0.0) new_a += tPi;

	return new_a;
}

void
//%----------------------------------------------------------------------------
//% Computes sinus of an interval
//%----------------------------------------------------------------------------
iSin(interval& w, interval& x)
{
	volatile real l = 0.0, u = 0.0;
	real inf, sup, xinf = x.inf(), xsup = x.sup();

	if (x.width() >= tPi) {
		hull(w, -1.0, 1.0);
	}
	else {
		l = NormalizeAngle(xinf);
		u = NormalizeAngle(xsup);

		if (u < l) {
			interval x1;
			interval w1, w2;

			x1 = interval{ l, tPi };

			iSin(w1, x1);

			x1 = interval{ 0.0, u };

			iSin(w2, x1);

			w = interval((w1.inf() < w2.inf()) ? w1.inf() : w2.inf(), (w1.sup() > w2.sup()) ? w1.sup() : w2.sup());
		}
		else {
			//% this is only valid for
			real a, b;
			RoundUp();
			a = sin(l);
			b = sin(u);
			RoundNear();

			if ((l <= hPi) && (hPi <= u))
				sup = 1.0;
			else
				sup = (a > b ? a : b);

			RoundUp();
			a = sin(l);
			b = sin(u);
			RoundNear();

			if ((l <= (3.0*hPi)) && ((3.0*hPi) <= u))
				inf = -1.0;
			else
				inf = (a < b ? a : b);
			w = interval{ inf, sup };
		}
	}
}

void
//%----------------------------------------------------------------------------
//% Scaling of x to lye in [0, 2Pi])
//%----------------------------------------------------------------------------
ScaleTo2Pi(interval& p)
{
	real q, x_inf, x_sup;
	if (p.inf() >= 0.0 && p.sup() <= hPi) { //% x in [0, pi/2]
		return;
	}
	q = p.inf() / tPi;
	q = floor(predR(q));
	q *= tPi;
	q = succR(q);
	if (q > 0.0) {
		q *= (1.0 + vnEps); //% vnEps
		q = succR(q);
	}
	x_inf = p.inf() - q;
	x_inf = predR(x_inf);
	if (x_inf < 0.0) {
		x_inf += tPi;
		x_inf = predR(x_inf);
	}
	q = p.sup() / tPi;
	q = floor(predR(q));
	q *= tPi;
	q = predR(q);
	if (q < 0.0) {
		q *= (1.0 + vnEps); //% vnEps
		q = predR(q);
	}
	x_sup = p.sup() - q;
	x_sup = succR(x_sup);
	if (x_sup > tPi) {
		x_sup -= tPi;
		x_sup = succR(x_sup);
	}
	p = interval{ x_inf, x_sup };
}

void
//%----------------------------------------------------------------------------
//% w = Sin(x)
//%----------------------------------------------------------------------------
iSin2(interval& w, interval& x)
{
	int q_inf, q_sup;
	real x_inf = x.inf(), x_sup = x.sup(), y_inf, y_sup;

	if (x.rad() >= tPi) {
		hull(w, -1.0, 1.0); //% = (-1) ^ (1); //HullRR(w, -1.0, 1.0);
		return;
	}
	//% Reduction to [0,2Pi]
	ScaleTo2Pi(x);
	if (x_sup >= x_inf + tPi) { //% security
		hull(w, -1.0, 1.0); //% w = (-1) ^ (1); //HullRR(w, -1.0, 1.0);
		return;
	}
	//% x is now in a range [0,2Pi]
	q_inf = Quadrant(x.inf());
	q_sup = Quadrant(x.sup());
	if ((q_inf == q_sup) && (x.sup() > x.inf() + Pi))
	{
		hull(w, -1.0, 1.0); //% , w = (-1) ^ (1); //HullRR(w, -1.0, 1.0);
		return;
	}
	switch ((q_sup << 2) + q_inf) {
	case 0:
	case 3:
	case 15:
		x_inf = sin(x.inf());
		x_sup = sin(x.sup());
		y_inf = rRoundDown(x_inf);
		y_sup = rRoundUp(x_sup);
		break;
	case 1:
	case 14:
		y_inf = -1.0;
		x_inf = sin(x_inf);
		x_sup = sin(x_sup);
		x_inf = rRoundUp(x_inf);
		x_sup = rRoundUp(x_sup);
		y_sup = ISLIB_MAX(x_inf, x_sup);
		break;
	case 2:
		y_inf = -1.0;
		x_sup = sin(x_sup);
		y_sup = rRoundUp(x_sup);
		break;
	case 4:
	case 11:
		y_sup = 1.0;
		x_inf = sin(x_inf);
		x_sup = sin(x_sup);
		x_inf = rRoundDown(x_inf);
		x_sup = rRoundDown(x_sup);
		y_inf = ISLIB_MIN(x_inf, x_sup);
		break;
	case 5:
	case 9:
	case 10:
		x_inf = sin(x_inf);
		x_sup = sin(x_sup);
		y_inf = rRoundDown(x_sup);
		y_sup = rRoundUp(x_inf);
		break;
	case 6:
	case 12:
		y_inf = -1.0; y_sup = 1.0; break;
	case 7:
		y_sup = 1.0;
		x_inf = sin(x_inf);
		y_inf = rRoundDown(x_inf);
		break;
	case 8:
		y_sup = 1.0;
		x_sup = sin(x_sup);
		y_inf = rRoundDown(x_sup);
		break;
	case 13:
		y_inf = -1.0;
		x_inf = sin(x_inf);
		y_sup = rRoundUp(x_inf);
		break;
	}
	if (y_inf < -1.0) y_inf = -1.0;     //%don't overestimate
	if (y_sup > 1.0) y_sup = 1.0;       //% dto. 
	hull(w, y_inf, y_sup); //% HullRR (w, y_inf, y_sup); 
}

void
//%----------------------------------------------------------------------------
//% w = Cos(x)
//%----------------------------------------------------------------------------
iCos(interval& w, interval& x)
{
	interval tmp;
	interval hPiIncl;

	interval{ hPi }.succ(hPiIncl);
	tmp = x + hPiIncl;
	iSin(w, tmp);
}

void
//%----------------------------------------------------------------------------
//% w = Tan(x) (tangens)
//%----------------------------------------------------------------------------
iTan(interval& w, interval& x)
{
	real x_inf = x.inf();
	real x_sup = x.sup();
	real y_inf, y_sup;
	int q_inf, q_sup;

	if (x.rad() >= tPi) {
		hull(w, NegInf, PosInf);
		return;
	}
	//% Reduction of x to [0, 2Pi]
	ScaleTo2Pi(x);
	if (x_sup >= x_inf + tPi) { //% security
		hull(w, NegInf, PosInf);
		return;
	}
	//% x_inf, x_sup are now in the range [0, 2Pi]
	//% Quadrants:
	//%  0 = [0,Pi/2]
	//%  1 = [Pi/2,Pi]
	//%  2 = [Pi,3Pi/2]
	//%  3 = [3Pi/2,2Pi]
	q_inf = Quadrant(x_inf);
	q_sup = Quadrant(x_sup);

	if ((q_inf == q_sup) && (x_sup > x_inf + Pi)) {
		hull(w, NegInf, PosInf);
		return;
	}
	switch ((q_sup << 2) + q_inf)
	{
	case 0:
	case 3:
	case 5:
	case 9:
	case 10:
	case 15:
		x_inf = tan(x_inf);
		x_sup = tan(x_sup);
		y_inf = rRoundDown(x_inf);
		y_sup = rRoundUp(x_sup);
		break;
	default:
		y_inf = NegInf;
		y_sup = PosInf;
		break;
	}
	hull(w, y_inf, y_sup);
}

void
//%----------------------------------------------------------------------------
//% w = Cot(x) (cotangens)
//%----------------------------------------------------------------------------
iCot(interval& w, interval& x)
{
	interval t1, t2;
	interval hPiIncl;
	interval(hPi).succ(hPiIncl);
	t1 = x + hPiIncl;
	iTan(t2, t1);
	w = -1.0 * t2;
}

void
//%----------------------------------------------------------------------------
//% w = Exp(x)
//%----------------------------------------------------------------------------
iExp(interval& w, interval& x)
{
	interval xx = x;
	real y_inf, y_sup;

	y_inf = std::nextafter(exp(xx.inf()), -ISLIB_INFINITY);
	y_sup = std::nextafter(exp(xx.sup()), ISLIB_INFINITY);
	if (y_inf < 0.0) y_inf = 0.0;
	hull(w, y_inf, y_sup);
}

void
//%----------------------------------------------------------------------------
//% w = Ln(x)
//%----------------------------------------------------------------------------
iLog(interval& w, interval& x)
{
	real y_inf, y_sup;

	if (x.inf() <= 0.0) {
		std::cout << "Log argument out of range"; /* BiasHullR (pR, & BiasNaN); */
	}
	else {
		y_inf = std::nextafter(log(x.inf()), -ISLIB_INFINITY);
		y_sup = std::nextafter(log(x.sup()), ISLIB_INFINITY);
		y_inf = rRoundDown(y_inf);
		y_sup = rRoundUp(y_sup);
		hull(w, y_inf, y_sup);
	}
}

void
//%----------------------------------------------------------------------------
//% w = Log10(x)
//%----------------------------------------------------------------------------
iLog10(interval& w, interval& x)
{
	interval t1;

	if (x.inf() <= 0.0) {
		std::cout << "Log10 argument out of range"; //% BiasHullR (pR, & BiasNaN);
	}
	else {
		iLog(t1, x);
		//% w = 
		//% BiasDivII (pR, & t1, & BiasLn10Incl);
	}
}

interval
//%----------------------------------------------------------------------------
//% w = Sqr(x)
//%----------------------------------------------------------------------------
iSqr(const interval& x)
{
	interval w;
	if (x.set_contains(0.0)) {
		real p1, p2;
		p1 = std::nextafter(x.inf()*x.inf(), ISLIB_INFINITY);
		p2 = std::nextafter(x.sup()*x.sup(), ISLIB_INFINITY);
		w = interval(0.0, ISLIB_MAX(p1, p2));
	}
	else {
		real p1, p2, p3, p4;
		p1 = std::nextafter(x.inf()*x.inf(), -ISLIB_INFINITY);
		p2 = std::nextafter(x.sup()*x.sup(), -ISLIB_INFINITY);
		p1 = std::nextafter(x.inf()*x.inf(), ISLIB_INFINITY);
		p2 = std::nextafter(x.sup()*x.sup(), ISLIB_INFINITY);
		w = interval(ISLIB_MIN(p1, p2), ISLIB_MAX(p3, p4));
	}
	return w;
}

interval
//%----------------------------------------------------------------------------
//% w = Sqrt(x)
//% Uses std::nextafter function to perform rounding
//%----------------------------------------------------------------------------
iSqrt(const interval& x)
{
	real infw, supw;
	interval w = interval::emptyset();
	if (x.inf() > 0.0) {
		infw = std::nextafter(sqrt(x.inf()), -ISLIB_INFINITY); //% RoundDown
		supw = std::nextafter(sqrt(x.sup()), ISLIB_INFINITY); //% RoundUp
		w = interval(infw, supw);
	}
	return w;
}