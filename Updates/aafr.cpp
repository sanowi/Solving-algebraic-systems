//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Class: Revised Affine Forms 
//%            Based on results of Kolev, Messine, Vu & Haroud, Hansen
//%
//%##########################################################################
//%

#include "../utils/stdafx.h"
#include "../interval/Round.h"
#include "multaux.h"
#include "aafr.h"

using namespace islib_interval;

int MULT = MULTS_BEST;// MULT_NEWBEST; //% selected multiplication

//%----------------------------------------------------------------------------
//% Constructor from a revised affine form from an affine form
//%----------------------------------------------------------------------------
AAFR::AAFR(const AAFR& a) : AAF(a) { _r = a._r; }

AAFR&
//%----------------------------------------------------------------------------
//% Copying operator
//%----------------------------------------------------------------------------
AAFR::operator = (const AAFR& a)
{
	if (this != &a) {
		_coeffs = a._coeffs;
		_indexes = a._indexes;
		_r = a._r;
	}
	return *this;
}

std::ostream&
//%----------------------------------------------------------------------------
//% Print a revised affine form
//%----------------------------------------------------------------------------
operator << (std::ostream& os, const AAFR& a)
{
	dvector	cfs = a._coeffs;
	cvector	idx = a._indexes;

	os.precision(12);
	os << cfs[0];
	for (int i = 1; i < cfs.size(); ++i) {
		if (cfs[i] > 0) {
			os << " + " << cfs[i];
		}
		else {
			os << " - " << ISLIB_ABS(cfs[i]);
		}
		os << "e_" << idx[i];
	}
	if (a._r != 0) {
		if (a._r > 0) {
			os << " + " << a._r << "r";
		}
		else {
			os << " - " << ISLIB_ABS(a._r) << "r";
		}
	}
	return os;
}

interval
//%----------------------------------------------------------------------------
//% Verified conversion from an affine form to an interval 
//%----------------------------------------------------------------------------
AAFR::reduce() const
{
	int ltemp = _indexes.size();
	real pr, inf, sup;

	pr = 0;
	RoundUp();
	for (int i = 1; i < ltemp; ++i) {
		pr += ISLIB_ABS(_coeffs[i]); //% radius of an affine form
	}
	inf = -(-_coeffs[0] + pr + _r);
	sup = _coeffs[0] + pr + _r;
	RoundNear();
	return interval{ inf, sup };
}

AAFR
//%----------------------------------------------------------------------------
//% Verfied addition of affine forms (Stolfi)
//%-----------------------------------------------------------------------------
operator + (const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size(); //% size of x
	int ny = y._indexes.size(); //% size of y
	real d = 0.0; //% rounding error
	cvector	xidx = x._indexes; //% indexes of x
	cvector	yidx = y._indexes, zidx; //% indexes of y
	dvector	xcfs = x._coeffs; //% coefficients of x
	dvector	ycfs = y._coeffs; //% coefficients of y

	zidx = AAF::idx_union(xidx, yidx); //% sum of sets of indexes
	ltemp = zidx.size(); //% cardinality of zidx

	dvector zcfs{ ltemp }; //% nzcfs - near mode
	dvector uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector dzcfs{ ltemp }; //% dzcfs - round-down mode

	//% Round near mode
	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = ycfs[iy++]; //% y.coeffs[ ib ] + 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++]; //% x.coeffs[ ia ] + 0
			continue;
		}
		zcfs[iz] = xcfs[ix++] + ycfs[iy++];
	}
	//% Round up mode
	RoundUp();
	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		dzcfs[iz] = -(-xcfs[ix] - ycfs[iy]);
		uzcfs[iz] = xcfs[ix++] + ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	d += x._r + y._r;
	//% Return to round near mode
	RoundNear();
	return AAFR{ zidx, zcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Verified addition of an affine form and a real number (on the right)
//%----------------------------------------------------------------------------
operator + (const AAFR& x, const real y)
{
	real zx, dx, ux, d = 0.0;
	cvector	xidx = x._indexes;
	dvector	xcfs = x._coeffs;

	//% Round near mode
	zx = xcfs[0] + y;
	//% Round up mode
	RoundUp();
	dx = -(-xcfs[0] - y);
	ux = xcfs[0] + y;
	d = x._r + ISLIB_MAX(ux - zx, zx - dx);
	//% Return to round near mode
	RoundNear();
	xcfs[0] = zx;
	return AAFR{ xidx, xcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Verified addition of a real number (on the left) and an affine form
//%----------------------------------------------------------------------------
operator + (const real y, const AAFR& x)
{
	real zx, dx, ux, d = 0.0;
	cvector	xidx = x._indexes;
	dvector	xcfs = x._coeffs;

	zx = xcfs[0] + y;
	RoundUp();
	dx = -(-xcfs[0] - y);
	ux = xcfs[0] + y;
	d = x._r + ISLIB_MAX(ux - zx, zx - dx);
	RoundNear();
	xcfs[0] = zx;
	return AAFR{ xidx, xcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Verified addition of an affine form and a real number (on the right)
//%----------------------------------------------------------------------------
operator + (const AAFR& x, const interval y)
{
	cvector	xidx = x._indexes;
	dvector	xcfs = x._coeffs;
	real d;

	RoundUp();
	xcfs[0] += y.mid(); //% changed 30.12.2018
	d = x._r + y.rad();
	RoundNear();
	return AAFR{ xidx, xcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Additive inverse (unary operator)
//%----------------------------------------------------------------------------
operator - (const AAFR& x)
{
	int n = x._indexes.size();
	real d = x._r;
	dvector	ncoeffs(n);

	for (int i = 0; i < n; ++i) {
		ncoeffs[i] = -(x._coeffs[i]);
	}
	return AAFR{ x._indexes, ncoeffs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Verified substraction of affine forms (Stolfi)
//%----------------------------------------------------------------------------
operator - (const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0.0; // d - rounding error
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx);
	ltemp = zidx.size();

	dvector zcfs{ ltemp }; //% nzcfs - near mode
	dvector uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector dzcfs{ ltemp }; //% dzcfs - round-down mode

	//% RoundNear mode
	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = -ycfs[iy++];  //% 0 - y.coeffs[ ib ]
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++];  //% x.coeffs[ ia ] - 0
			continue;
		}
		zcfs[iz] = xcfs[ix++] - ycfs[iy++];
	}
	RoundUp();
	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		dzcfs[iz] = -(-xcfs[ix] + ycfs[iy]);
		uzcfs[iz] = xcfs[ix++] - ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	d += x._r + y._r;
	RoundNear();
	return AAFR{ zidx, zcfs, d };
}

AAFR
//%-----------------------------------------------------------------------------
//% Verified subtraction of a real number from an affine form
//%----------------------------------------------------------------------------
operator - (const AAFR& x, const real y)
{
	real zx, dx, ux, d = 0.0;
	cvector	xidx = x._indexes;
	dvector	xcfs = x._coeffs;

	zx = xcfs[0] - y;
	RoundUp();
	dx = -(-xcfs[0] + y);
	ux = xcfs[0] - y;
	d = x._r + ISLIB_MAX(ux - zx, zx - dx);
	RoundNear();
	xcfs[0] = zx;
	return AAFR{ xidx, xcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Verified subtraction of an affine form from a real number
//%----------------------------------------------------------------------------
operator - (const real y, const AAFR& x)
{
	real zx, dx, ux;
	real d = 0.0;
	cvector	xidx = x._indexes;
	dvector	xcfs = x._coeffs;
	int ltemp = xidx.size();

	for (int i = 1; i < ltemp; ++i) {
		xcfs[i] = -xcfs[i];
	}
	zx = y - xcfs[0];
	RoundUp();
	dx = -(-y + xcfs[0]);
	ux = y - xcfs[0];
	d = x._r + ISLIB_MAX(ux - zx, zx - dx);
	xcfs[0] = zx;
	RoundNear();
	return AAFR{ xidx, xcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Multiplication of affine forms based on Miyaima's result from the paper
//% "On the best multiplication of affine forms"
//% (minimum-error Chebyshev approximation)
//%----------------------------------------------------------------------------
AAFR::bestMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	interval d;
	interval rx, ry;
	real ax = 0, ay = 0, sx, sy, ez;
	real dmin, dmax;
	cvector xidx = x._indexes;
	cvector yidx = y._indexes, zidx;
	dvector xcfs = x._coeffs;
	dvector ycfs = y._coeffs;

	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) {
		return xcfs[0] * ycfs[0];
	}
	if (rx.rad() == 0.0) {
		return xcfs[0] * y;
	}
	if (ry.rad() == 0.0) {
		return x * ycfs[0];
	}

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp };

	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}

	dMinMax(0.0, 0.0, xcfs2, ycfs2, dmin, dmax);

	d = interval{ dmin, dmax };
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	sx = ISLIB_ABS(xcfs[0]);
	sy = ISLIB_ABS(ycfs[0]);
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
		sx += ISLIB_ABS(xcfs2[i]);
		sy += ISLIB_ABS(ycfs2[i]);
	}
	ez = d.rad() + x._r * sy + y._r * sx + x._r * y._r;
	return AAFR{ zidx, zcfs, ez };
} //% >>> bestMult <<<

AAFR
//%----------------------------------------------------------------------------
//% O(n log(n)) algorithm for Chebyshev multiplication of revised affine forms
//% This multiplication is verified to some extent.
//%----------------------------------------------------------------------------
AAFR::newBestMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	interval d, rx, ry;
	real ax = 0, ay = 0, ez, dmin, dmax;
	cvector	xidx = x._indexes, yidx = y._indexes, zidx;
	dvector xcfs = x._coeffs, ycfs = y._coeffs;

	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) {
		return xcfs[0] * ycfs[0];
	}
	if (rx.rad() == 0.0) {
		return xcfs[0] * y;
	}
	if (ry.rad() == 0.0) {
		return x * ycfs[0];
	}

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp };

	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}

	int lltemp = ltemp - 1;
	if (x._r != 0.0) lltemp++;
	if (y._r != 0.0) lltemp++;
	dvector xcfs3{ lltemp }, ycfs3{ lltemp };
	for (int i = 1; i < ltemp; ++i) {
		xcfs3[i - 1] = xcfs2[i];
		ycfs3[i - 1] = ycfs2[i];
	}
	int ii = ltemp - 1;
	if (x._r != 0.0) {
		xcfs3[ii] = x._r;  ycfs3[ii] = 0.0; ++ii;
	}
	if (y._r != 0.0) {
		xcfs3[ii] = 0.0; ycfs3[ii] = y._r;
	}
	//% Compute minimum and maximum of f(x,y)=x*y on the joint range
	dmm(xcfs3, ycfs3, dmin, dmax);
	d = interval{ dmin, dmax }; //% Computing with round to near
	
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	RoundUp();
	dzcfs[0] = -(-(xcfs[0] * ycfs[0]) - d.mid());
	for (int i = 1; i < ltemp; ++i) {
		//% By reversing the sign the rounding down is achieved
		dzcfs[i] = -(-(xcfs[0] * ycfs2[i]) - (ycfs[0] * xcfs2[i]));
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	real dd = 0.0;
	for (int i = 0; i < ltemp; ++i) {
		dd += ISLIB_MAX(ISLIB_ABS(zcfs[i] - dzcfs[i]), ISLIB_ABS(uzcfs[i] - zcfs[i]));
	}
	ez = d.rad() + x._r * ISLIB_ABS(ycfs[0]) + y._r * ISLIB_ABS(xcfs[0]) + dd;
	RoundNear();
	return AAFR{ zidx, zcfs, ez };
} //% >>> newBestMult <<<

AAFR
//%----------------------------------------------------------------------------
//% O(n) algorithm for Chebyshev multiplication of revised affine forms.
//% The asymptotic time comlexity of the algorithm is O(n), but it generally
//% is slower than its O(nlog(n)) version.
//%----------------------------------------------------------------------------
AAFR::newBestMultON(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	interval d, rx, ry;
	real ax = 0, ay = 0, ez, dmin, dmax;
	cvector	xidx = x._indexes, yidx = y._indexes, zidx;
	dvector	 xcfs = x._coeffs, ycfs = y._coeffs;

	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) {
		return xcfs[0] * ycfs[0];
	}
	if (rx.rad() == 0.0) {
		return xcfs[0] * y;
	}
	if (ry.rad() == 0.0) {
		return x * ycfs[0];
	}

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp };

	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	int lltemp = ltemp - 1;
	if (x._r != 0.0) lltemp++;
	if (y._r != 0.0) lltemp++;
	dvector xcfs3{ lltemp }, ycfs3{ lltemp };
	for (int i = 1; i < ltemp; ++i) {
		xcfs3[i - 1] = xcfs2[i];
		ycfs3[i - 1] = ycfs2[i];
	}
	int ii = ltemp - 1;
	if (x._r != 0.0) {
		xcfs3[ii] = x._r;  ycfs3[ii] = 0.0; ++ii;
	}
	if (y._r != 0.0) {
		xcfs3[ii] = 0.0; ycfs3[ii] = y._r;
	}
	int l;
	real vx0, vy0;
	dvector hx, hy, hs;
	dmm_slopes(xcfs3, ycfs3, l, vx0, vy0, hx, hy, hs);
	dmm_test_max(xcfs3, ycfs3, dmax, l, vx0, vy0, hx, hy, hs);
	dmm_test_min(xcfs3, ycfs3, dmin, l, vx0, vy0, hx, hy, hs);
	d = interval{ dmin, dmax };
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	RoundUp();
	dzcfs[0] = -(-(xcfs[0] * ycfs[0]) - d.mid());
	for (int i = 1; i < ltemp; ++i) {
		//% By reversing the sign the rounding down is achieved
		dzcfs[i] = -(-(xcfs[0] * ycfs2[i]) - (ycfs[0] * xcfs2[i]));
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	real dd = 0.0;
	for (int i = 0; i < ltemp; ++i) {
		dd += ISLIB_MAX(ISLIB_ABS(zcfs[i] - dzcfs[i]), ISLIB_ABS(uzcfs[i] - zcfs[i]));
	}
	ez = d.rad() + x._r * ISLIB_ABS(ycfs[0]) + y._r * ISLIB_ABS(xcfs[0]) + dd;
	RoundNear();
	return AAFR{ zidx, zcfs, ez };
} //% >>> newBestMultON <<<

AAFR
//%----------------------------------------------------------------------------
//% Verified standard multiplication of affine forms Kolev's version
//%----------------------------------------------------------------------------
AAFR::stdMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0, s = 0, s2 = 0, s3 = 0;
	real sup = 0, sdown = 0, sabs = 0, es, ds, g;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of sets of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; //% nzcfs - near mode
	dvector	uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; //% dzcfs - round-down mode

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * s; //% x[ 0 ] * y[ 0 ] + 0.5 * sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		sabs += ISLIB_ABS(g); //% suma( | x[ i ] * y[ i ] | ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sdown; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	es = 0.5 * sabs;
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); // v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); // u
		ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sup; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	ds = ISLIB_ABS(x._coeffs[0]) * y._r + ISLIB_ABS(y._coeffs[0]) * x._r + (s2 + x._r) * (s3 + y._r) + d;
	RoundNear();
	ds -= es;
	return AAFR{ zidx, zcfs, ds };
} //% >> stdMult <<<

AAFR
//%----------------------------------------------------------------------------
//% Verified trivial multiplication of affine forms Kolev's version
//%----------------------------------------------------------------------------
AAFR::quickMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0, s = 0, s2 = 0, s3 = 0;
	real sup = 0, sdown = 0, sabs = 0, ds, g;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of sets of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; //% nzcfs - near mode
	dvector	uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; //% dzcfs - round-down mode

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		sabs += ISLIB_ABS(g); //%  suma( | x[ i ] * y[ i ] | ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); //% u
		ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	ds = x._r * y._r + y._r * (ISLIB_ABS(x._coeffs[0]) + s2) + x._r * (ISLIB_ABS(y._coeffs[0]) + s3) + s2 * s3 + d;
	RoundNear(); //%!!!
	return AAFR{ zidx, zcfs, ds };
} //% >>> quickMult <<<

AAFR
//%----------------------------------------------------------------------------
//% Modified standard multiplication of affine forms by Miyajima
//%----------------------------------------------------------------------------
AAFR::miyMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0.0, s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
	real sup = 0.0, sdown = 0.0, sabs = 0.0, es, ds, g;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;
	
	zidx = AAF::idx_union(xidx, yidx); // sum of sets of indexes
	ltemp = zidx.size();

	dvector	zcfs{ ltemp }; // zcfs - near mode
	dvector	uzcfs{ ltemp }; // uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; // dzcfs - round-down mode
	
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * s; //% x[ 0 ] * y[ 0 ] + 0.5 * sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sdown; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sup += g;
		sabs += ISLIB_ABS(g); //% suma( | x[ i ] * y[ i ] | ) - round up
		s3 += ISLIB_ABS(ycfs[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); //% u
		ix++;
	}
	es = 0.5 * sabs;
	uzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sup; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	dvector xcfs2{ ltemp }, ycfs2{ ltemp };
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	//% computing sum_{i<j}|x_iy_j+x_jy_i|
	for (int i = 1; i < ltemp; ++i) {
		for (int j = 1; j < i; ++j) {
			s4 += ISLIB_ABS(xcfs2[i] * ycfs2[j] + xcfs2[j] * ycfs2[i]);
		}
	}
	ds = x._r * y._r + y._r * (ISLIB_ABS(x._coeffs[0]) + s2) + x._r * (ISLIB_ABS(y._coeffs[0]) + s3) +
		s4 + es + d;
	RoundNear(); //%!!!
	return AAFR{ zidx, zcfs, ds };
} //% >>> miyMult <<<

AAFR
//%----------------------------------------------------------------------------
//% Multiplication based on tangent plane
//%----------------------------------------------------------------------------
AAFR::AASEMult2(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp, ltemp2, nn;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	int nx2 = nx, ny2 = ny;
	real d = 0, s = 0, s2 = 0, s3 = 0;
	real sup = 0, sdown = 0, sabs = 0, es, g;
	cvector	xidx = x._indexes, xidx2, yidx2;
	cvector	yidx = y._indexes, zidx, zidx2;
	dvector	xcfs = x._coeffs, xcfs2;
	dvector	ycfs = y._coeffs, ycfs2, zcfs2;

	xidx2 = xidx; xcfs2 = xcfs;
	yidx2 = yidx; ycfs2 = ycfs;
	zidx = AAFR::idx_union(xidx, yidx);
	ltemp = zidx.size();
	dvector zcfs{ ltemp };

	//% if ncessary (i.e. accumulative errora are non-zero), we expand the revised affine form
	nn = AAFR::getlast();
	if (x.r() != 0.0) {
		nx2++;
		xidx2 = cvector{ nx2 };
		xcfs2 = dvector{ nx2 };
		for (int i = 0; i < nx; ++i) {
			xidx2[i] = xidx[i];
			xcfs2[i] = xcfs[i];
		}
		nn++;
		xidx2[nx] = nn;
		xidx2[nx] = (int)x.r();
	}
	if (y.r() != 0.0) {
		ny2++;
		yidx2 = cvector{ ny2 };
		ycfs2 = dvector{ ny2 };
		for (int i = 0; i < ny; ++i) {
			yidx2[i] = yidx[i];
			ycfs2[i] = ycfs[i];
		}
		nn++;
		yidx2[ny] = nn;
		yidx2[ny] = (int)y.r();
	}
	zidx2 = AAF::idx_union(xidx2, yidx2); //% sum of sets of indexes
	ltemp2 = zidx2.size();
	zcfs2 = dvector{ ltemp2 };
	zcfs[0] = xidx[0] * yidx[0];
	for (ix = 1, iy = 1, iz = 1; iz < ltemp2; ++iz) {
		if (ix == nx2 || xidx2[ix] != zidx2[iz]) {
			zcfs2[iz] = xcfs2[0] * ycfs2[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny2 || yidx2[iy] != zidx2[iz]) {
			zcfs2[iz] = xcfs2[ix++] * ycfs2[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs2[iz] = xcfs2[ix++] * ycfs2[0] + xcfs2[0] * ycfs2[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	sdown = 0.0; sabs = 0.0;
	for (ix = 1, iy = 1, iz = 1; iz < ltemp2; ++iz) {
		if (ix == nx2 || xidx2[ix] != zidx2[iz]) {
			iy++;
			continue;
		}
		if (iy == ny2 || yidx2[iy] != zidx2[iz]) {
			ix++;
			continue;
		}
		g = xcfs2[ix] * ycfs2[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		sabs += ISLIB_ABS(g); //% suma( |x[ i ] * y[ i ]| )
		ix++;
		iy++;
	}
	zcfs[0] += sdown;
	/*for (ix = 1, iz = 1; iz < ltemp2; ++iz) {
		if (ix == nx2 || xidx2[ix] != zidx2[iz]) {
			continue;
		}
		for (iy = 1, iz1 = 1; iz1 <ltemp2; iz1++) {
			if (iy == ny2 || yidx2[iy] != zidx2[iz1]) {
				continue;
			}
			if (ix != iy) {
				s += ISLIB_ABS(xcfs2[ix] * ycfs2[iy]);
			}
			iy++;
		}
		ix++;
	}*/
	for (ix = 1, iy = 1, iz = 1; iz < ltemp2; ++iz) {
		if (ix == nx2 || xidx2[ix] != zidx2[iz]) {
			s3 += ISLIB_ABS(ycfs2[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny2 || yidx2[iy] != zidx2[iz]) {
			s2 += ISLIB_ABS(xcfs2[ix]); //% u - round up
			ix++; continue;
		}
		sup += xcfs2[ix] * ycfs2[iy];
		s3 += ISLIB_ABS(ycfs2[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs2[ix]); //% u
		ix++;
	}
	es = s3 * s2 - sabs;
	return AAFR{ zidx, zcfs, es };
}

AAFR
//%----------------------------------------------------------------------------
//% Based on the paper:
//% A refined affine approximation method of multiplication for range analysis 
//% in word-length optimization
//% AAFR operator * (const AAFR& x, const AAFR& y)
//% This multiplication is QUESTIONABLE!!!
//%----------------------------------------------------------------------------
AAFR::AASEMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0, s = 0, s2 = 0, s3 = 0;
	real sup = 0, sdown = 0, sabs = 0, es, ds, g;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of sets of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; //% nzcfs - near mode
	dvector	uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; //% dzcfs - round-down mode

	if (ltemp == 2) {
		return AAFR::newBestMult(x, y);
	}
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0] + s; //% x[ 0 ] * y[ 0 ] + sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		sabs += ISLIB_ABS(g); //% suma( | x[ i ] * y[ i ] | ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + sdown; //% x[ 0 ] * y[ 0 ] + suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	es = sabs;
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); //% u
		ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + sup; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	ds = ISLIB_ABS(x._coeffs[0]) * y._r + ISLIB_ABS(y._coeffs[0]) * x._r + (s2 + x._r) * (s3 + y._r) + d;
	RoundNear(); //%!!!
	ds -= es;
	return AAFR{ zidx, zcfs, ds };
}

AAFR
//%----------------------------------------------------------------------------
//% Modified standard multiplication of affine forms Rump - Kashiwagi version
//%----------------------------------------------------------------------------
AAFR::rukasMult(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0.0, s = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
	real sup = 0.0, sdown = 0.0, seplus = 0.0, seminus = 0.0, ds, g;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); // sum of sets of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; // nzcfs - near mode
	dvector	uzcfs{ ltemp }; // uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; // dzcfs - round-down mode

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * s; //% x[ 0 ] * y[ 0 ] + 0.5 * sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  // 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		if (g < 0.0) {
			seminus += g; //% suma( x[ i ] * y[ i ] < 0) - round down
		}
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sdown; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		g = xcfs[ix] * ycfs[iy];
		if (g >= 0.0) {
			seplus += g; //% suma( x[ i ] * y[ i ] >= 0) - round down
		}
		sup += g;
		s3 += ISLIB_ABS(ycfs[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); //% u
		ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sup; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}

	dvector xcfs2{ ltemp }, ycfs2{ ltemp };
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	for (int i = 1; i < ltemp; ++i) {
		for (int j = 1; j < i; ++j) {
			s4 += ISLIB_ABS(xcfs2[i] * ycfs2[j] + xcfs2[j] * ycfs2[i]);
		}
	}
	ds = x._r * y._r + y._r * (ISLIB_ABS(x._coeffs[0]) + s2) + x._r * (ISLIB_ABS(y._coeffs[0]) + s3) +
		0.5 * ISLIB_MAX(seplus, -seminus) + s4 + d;
	RoundNear(); //%!!!
	return AAFR{ zidx, zcfs, ds };
} //% >>> rukasMult <<<

AAFR
//%-----------------------------------------------------------------------------
//% Trivial multiplication of affine forms (Kolev's version)
//%----------------------------------------------------------------------------
AAFR::quickMultKolev(const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	real d = 0, s = 0, s2 = 0, s3 = 0;
	real sup = 0, sdown = 0, ds, g;
	real xnorm, ynorm, ss1, ss2;
	cvector	xidx = x._indexes;
	cvector	yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs;
	dvector	ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of sets of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; //% nzcfs - near mode
	dvector	uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector	dzcfs{ ltemp }; //% dzcfs - round-down mode

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum( x_i * y_i ) - round near
	}
	zcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * sum( x[ i ] * y[ i ] ) - round near
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0  - round near
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% 0 + x[ i ] * y[ 0 ] - round near
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ 0 ] * y[ i ] + x[ i ] * y[ 0 ] - round near
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		g = xcfs[ix] * ycfs[iy];
		sdown += g; //% suma( x[ i ] * y[ i ] ) - round down
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round down
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] - round down
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ia]*0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++]; //% x[ i ]*y[0] + x[ i ]*y[0] - round down
	}
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]); //% v - round up
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]); //% u - round up
			ix++; continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); //% v
		iy++;
		s2 += ISLIB_ABS(xcfs[ix]); //% u
		ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0]; //% x[ 0 ] * y[ 0 ] + 0.5 * suma( x[ i ] * y[ i ] ) - round up mode
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% x[ 0 ] * y[ i ] + 0 - round up mode
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			uzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
			continue;
		}
		uzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	dvector xcfs2{ ltemp }, ycfs2{ ltemp };
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exists
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { //% y index does not exists
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	xnorm = norm(xcfs2);
	ynorm = norm(ycfs2);
	if (xnorm != 0) {
		xcfs2 = xcfs2 * (1.0 / xnorm);
	}
	if (ynorm != 0) {
		ycfs2 = ycfs2 * (1.0 / ynorm);
	}
	ss1 = 0.0;
	ss2 = 0.0;
	for (int i = 0; i < ltemp; ++i) {
		ss1 += ISLIB_ABS(xcfs2[i] + ycfs2[i]);
		ss2 += ISLIB_ABS(xcfs2[i] - ycfs2[i]);
	}
	ds = x._r * y._r + y._r * (ISLIB_ABS(x._coeffs[0]) + s2) + x._r * (ISLIB_ABS(y._coeffs[0]) + s3) +
		0.25*xnorm*ynorm*ISLIB_MAX(ss1*ss1, ss2*ss2) + d;
	RoundNear(); //%!!!
	return AAFR{ zidx, zcfs, ds };
} //% >>> quickMultKolev <<<


AAFR
//%----------------------------------------------------------------------------
//% Multiplication of revised affine forms
//%----------------------------------------------------------------------------
operator * (const AAFR& x, const AAFR& y)
{
	AAFR z;
	switch (MULT) {
	case MULT_BEST:
		z = AAFR::bestMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD0!!!" << std::endl;
		return z;
	case MULT_NEWBEST:
		z = AAFR::newBestMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD1!!!" << std::endl;
		return z;
	case MULT_NEWBESTON:
		z = AAFR::newBestMultON(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD2!!!" << std::endl;
		return z;
	case MULT_STD:
		z = AAFR::stdMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD3!!!" << std::endl;
		return z;
	case MULT_QUICK:
		z = AAFR::quickMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD4!!!" << std::endl;
		return z;
	case MULT_MIY:
		z = AAFR::miyMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD5!!!" << std::endl;
		return z;
	case MULT_RUKAS:
		z = AAFR::rukasMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD6!!!" << std::endl;
		return z;
	case MULTS_AASEE:
		z = AAFR::AASEMult(x, y);
		if (z.reduce().is_empty()) std::cout << "BLAD7!!!" << std::endl;
		return z;
	default:
		assert(false);
		return 0;
	}
}

AAFR
//%----------------------------------------------------------------------------
//% Multiplication of a revised affine form and a real number (on the right)
//%----------------------------------------------------------------------------
operator * (const real x, const AAFR& y)
{
	int ltemp = y._indexes.size();
	real d = 0.0;
	cvector zidx = y._indexes;
	dvector zcfs = y._coeffs;
	dvector zx{ ltemp }, dz{ ltemp }, uz{ ltemp };

	for (int i = 0; i < ltemp; ++i) {
		zx[i] = 0;
		if (zcfs[i] != 0 && x != 0.0) {
			zx[i] = x * zcfs[i];
		}
	}
	RoundDown();
	for (int i = 0; i < ltemp; ++i) {
		dz[i] = 0.0;
		if (zcfs[i] != 0 && x != 0.0) {
			dz[i] = x * zcfs[i];
		}
	}
	RoundUp();
	for (int i = 0; i < ltemp; ++i) {
		uz[i] = 0;
		if (zcfs[i] != 0 && x != 0) {
			uz[i] = x * zcfs[i];
		}
		d += ISLIB_MAX(uz[i] - zx[i], zx[i] - dz[i]);
		zcfs[i] = zx[i];
	}
	d += ISLIB_ABS(x) * y._r;
	RoundNear();
	return AAFR{ zidx, zcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Multiplication of a revised affine form and a real number (on the left)
//%----------------------------------------------------------------------------
operator * (const AAFR& y, const real x)
{
	int ltemp = y._indexes.size();
	real d = 0.0;
	cvector zidx = y._indexes;
	dvector zcfs = y._coeffs;
	dvector zx{ ltemp }, dz{ ltemp }, uz{ ltemp };

	for (int i = 0; i < ltemp; ++i) {
		zx[i] = 0;
		if (zcfs[i] != 0 && x != 0.0) {
			zx[i] = zcfs[i] * x;
		}
	}
	RoundDown();
	for (int i = 0; i < ltemp; ++i) {
		dz[i] = 0.0;
		if (zcfs[i] != 0 && x != 0.0) {
			dz[i] = zcfs[i] * x;
		}
	}
	RoundUp();
	for (int i = 0; i < ltemp; ++i) {
		uz[i] = 0;
		if (zcfs[i] != 0 && x != 0) {
			uz[i] = zcfs[i] * x;
		}
		d += ISLIB_MAX(uz[i] - zx[i], zx[i] - dz[i]);
		zcfs[i] = zx[i];
	}
	d += ISLIB_ABS(x) * y._r;
	RoundNear();
	return AAFR{ zidx, zcfs, d };
}

AAFR
//%----------------------------------------------------------------------------
//% Trivial multiplication of revised affine forms
//%----------------------------------------------------------------------------
AAFR::quickmult(const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = _indexes.size();
	int ny = y._indexes.size();
	real d, sx, sy;
	cvector xidx = _indexes;
	cvector yidx = y._indexes, zidx;
	dvector xcfs = _coeffs;
	dvector ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();
	dvector zcfs(ltemp); //% nzcfs - near mode
	dvector xcfs2{ ltemp }, ycfs2{ ltemp };
	dvector einf{ ltemp }, esup{ ltemp }; //% dummy variables for min and max

	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			ycfs2[iz] = ycfs[iy++]; //% x index does not exists
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			xcfs2[iz] = xcfs[ix++]; //% y index does not exists
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * _r * y._r;
	sx = 0;
	sy = 0;
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
		sx += ISLIB_ABS(xcfs2[i]);
		sy += ISLIB_ABS(ycfs2[i]);
	}
	d = sx * sy + _r * (ISLIB_ABS(ycfs[0]) + sy) + y._r * (ISLIB_ABS(xcfs[0]) + sx) + 0.5 * _r * y._r;
	return AAFR{ zidx, zcfs, d };
} //% >>> quickmult <<<

bool
//%----------------------------------------------------------------------------
//% Equivalence of an affine form and a real number
//%----------------------------------------------------------------------------
operator == (const AAFR& a, const real b)
{
	if (a.reduce() == b) return true;
	return false;
}

AAFR
//%----------------------------------------------------------------------------
//% Rreciprocal of an affine form. Version: Kolev, Miyajima approximation 
//% using linear function: f(x) = ax + b + d*eps_new which is an "average"  
//% of tangent and secant lines.
//%----------------------------------------------------------------------------
AAFR::reciprocal() const
{
	int i, ltemp = _indexes.size();
	interval c = reduce();
	real cinf, csup, ab, a, x, x0, b1, b2, b, rw;
	cvector ridx = _indexes;
	dvector rcfs{ ltemp };

	RoundUp();
	cinf = c.inf(); //% c is as range of reduced affine form
	csup = c.sup();
	ab = cinf * csup;
	a = -1.0 / ab; //% slope (gradient)
	x = sqrt(ab);
	x0 = (c > 0.0) ? x : -x; //% x0 - tangent point
	b1 = 1.0 / cinf - a * cinf; //% b1 - secant intercept
	RoundDown();
	b2 = 1.0 / x0 - a * x0; //% b2 - tangent intercept
	RoundUp();
	b = b2 + 0.5 * (b1 - b2);  //% b = 0.5 * ( b2 + b1 ) - mid( [b1, b2] )
	rw = ISLIB_ABS(a) * _r + ISLIB_ABS(b - b2);  //% rw = b - b2 = 0.5*( b1 - b2 ) - rad( [ b1, b2 ] ) //b1 - b;
	RoundNear();
	rcfs[0] = a * _coeffs[0] + b;
	for (i = 1; i < ltemp; ++i) {
		rcfs[i] = a * _coeffs[i];
	}
	return AAFR{ ridx, rcfs, rw };
}

AAFR
//%----------------------------------------------------------------------------
//% Rreciprocal of a revised affine form (min-range approximation)
//%----------------------------------------------------------------------------
AAFR::reciprocal_mr() const
{
	int i, ltemp = _indexes.size();
	interval c = reduce();
	real cinf, csup, ab, a, b1, b2, b, rw;
	cvector ridx = _indexes;
	dvector rcfs{ ltemp };

	RoundUp();
	cinf = c.inf(); //% c is as range of reduced affine form
	csup = c.sup();
	ab = cinf * csup;
	a = -1.0 / (csup*csup); //% slope (gradient)
	b1 = (cinf*cinf + csup * csup) / (cinf*csup*csup); //% b1 - secant intercept
	RoundDown();
	b2 = 2.0 / csup; //% b2 - tangent intercept
	RoundUp();
	b = b2 + 0.5 * (b1 - b2);  //% b = 0.5 * ( b2 + b1 ) - mid( [b1, b2] )
	rw = ISLIB_ABS(a) * _r + ISLIB_ABS(b - b2);  //% rw = b - b2 = 0.5*( b1 - b2 ) - rad( [ b1, b2 ] ) //b1 - b;
	RoundNear();
	rcfs[0] = a * _coeffs[0] + b;
	for (i = 1; i < ltemp; ++i) {
		rcfs[i] = a * _coeffs[i];
	}
	return AAFR{ ridx, rcfs, rw };
}

AAFR
//%----------------------------------------------------------------------------
//% Reciprocal of a revised affine form. Rounding errors are not taken into
//% account.
//%----------------------------------------------------------------------------
AAFR::reciprocal2()
{
	int i, ltemp = _indexes.size();
	interval c = reduce();
	real cinf, csup, ab, rw;
	cvector ridx = _indexes;
	dvector rcfs{ ltemp };

	//% RoundUp();
	cinf = c.inf();
	csup = c.sup();
	ab = cinf * csup;
	for (i = 0; i < ltemp; ++i) {
		rcfs[i] = _coeffs[i] / ab;
	}
	rw = 0; //% for now
	RoundNear();
	return AAFR{ ridx, rcfs, rw };
}

AAFR
//%----------------------------------------------------------------------------
//% Division of revised affine forms. Version: Kolev, El Owny (see El-Owny).
//% Rounding errors are not taken into account.
//%----------------------------------------------------------------------------
operator / (const AAFR& x, const AAFR& y)
{
	int ix, iy, iz, ltemp;
	int nx = x._indexes.size();
	int ny = y._indexes.size();
	interval xi, yi;
	real c, rp;
	cvector xidx = x._indexes;
	cvector yidx = y._indexes, zidx, pidx;
	dvector xcfs = x._coeffs;
	dvector ycfs = y._coeffs, pcfs;
	AAFR p, ry, pry;

	pidx = AAF::idx_union(xidx, yidx);
	ltemp = pidx.size();
	pcfs = dvector{ ltemp };
	RoundUp();
	xi = x.reduce();
	yi = y.reduce();
	if (xi.rad() == 0.0) {
		return xi.mid() * y.reciprocal();
	}
	if (yi.rad() == 0.0) {
		return (1.0 / yi.mid() * x);
	}
	c = x._coeffs[0] / y._coeffs[0];
	pcfs[0] = 0;
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != pidx[iz]) {
			pcfs[iz] = -c * ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != pidx[iz]) {
			pcfs[iz] = xcfs[ix++];
			continue;
		}
		pcfs[iz] = xcfs[ix++] - c * ycfs[iy++];
	}
	rp = x._r + ISLIB_ABS(c) * y._r;
	p = AAFR{ pidx, pcfs, rp };
	ry = y.reciprocal();
	pry = p * ry;
	pry._coeffs[0] += c;
	RoundNear();
	return pry;
}

void
//%----------------------------------------------------------------------------
//% Incrementation by a revised affine form
//%----------------------------------------------------------------------------
operator += (AAFR& x, const AAFR& y)
{
	x = x + y;
}

void
//%----------------------------------------------------------------------------
//% Incrementation by a real number
//%----------------------------------------------------------------------------
operator += (AAFR& x, const real y)
{
	x = x + y;
}

void
//%----------------------------------------------------------------------------
//% Blow up a revised affine form
//%----------------------------------------------------------------------------
AAFR::blow(real eps, real mi)
{
	interval x = reduce();
	real r = x.rad();
	RoundUp();
	if (r > 0) {
		_r += r * eps;
	}
	else {
		_r += interval{ predR(x.inf()), succR(x.sup()) }.rad();
	}
	_r += mi;
	RoundNear();
}

AAFR
//%----------------------------------------------------------------------------
//% Square root of a revised affine form
//%----------------------------------------------------------------------------
AAFR::sqrtf()
{
	int i, ltemp = _indexes.size();
	real a, b, alpha, beta, delta, db, dmin, dmax;
	real fad, fau, fbd, fbu, aad, aau, abd, abu, gau, dfau;
	interval xi = reduce();
	cvector zidx = _indexes;
	dvector zcfs{ ltemp };

	a = xi.inf();
	b = xi.sup();
	if (xi.rad() == 0) {
		return AAFR{ sqrt(xi.mid()) };
	}
	RoundUp(); //% rounding up
	fau = sqrt(a);
	fbu = sqrt(b);
	alpha = (fbu - fau) / (b - a);
	aau = alpha * a;
	abu = alpha * b;
	gau = 1.0 / (4.0 * alpha);
	dfau = 1.0 / (2.0 * sqrt(a));

	RoundDown();  //% rounding down
	fad = sqrt(a);
	fbd = sqrt(b);
	aad = alpha * a;
	abd = alpha * b;

	RoundNear();  //% rounding to near
	db = fbd - abu;
	if (alpha > dfau) {
		dmin = db;
		dmax = fau - abd;
	}
	else {
		dmin = ISLIB_MIN(fad - aau, db);
		dmax = gau;
	}
	beta = interval{ dmin, dmax }.mid();
	delta = ISLIB_ABS(alpha)*_r + interval{ dmin, dmax }.rad();
	for (i = 0; i < ltemp; ++i) {
		zcfs[i] = alpha * _coeffs[i];
	}
	zcfs[0] += beta;
	return AAFR{ zidx, zcfs, delta };
}

AAFR
//%----------------------------------------------------------------------------
//% Square of a revised affine form (Chebyshev approximation)
//%----------------------------------------------------------------------------
AAFR::sqr()
{
	int ltemp = _indexes.size();
	real cinf, csup, b1, b2, rz, a;
	interval c;
	cvector zidx = _indexes;
	dvector zcfs{ ltemp };

	c = reduce();
	cinf = c.inf();
	csup = c.sup();
	a = cinf + csup;
	b1 = -(cinf * csup);
	b2 = -(a * a) / 4.0;
	for (int i = 0; i < ltemp; ++i) {
		zcfs[i] = _coeffs[i] * a;
	}
	zcfs[0] += (b1 + b2) / 2.0;
	rz = ISLIB_ABS(a) * _r + ISLIB_ABS(b1 - b2) / 2.0;
	return AAFR{ zidx, zcfs, rz };
}

AAFR
//%=----------------------------------------------------------------------------
//% Cube of a revised affine form 
//%----------------------------------------------------------------------------
AAFR::cube()
{
	int ltemp = _indexes.size();
	real a, b, b1, b2, d, x0, rz;
	interval c;
	dvector zcfs{ ltemp };
	cvector zidx = _indexes;

	c = reduce();
	a = c.inf();
	b = c.sup();
	d = (a + b) * (a + b) - a * b;
	x0 = 0.0;
	if (a > 0.0) x0 = sqrt(d / 3.0);
	if (b < 0.0) x0 = -sqrt(d / 3.0);
	b1 = -a * b * (a + b);
	b2 = -d * x0 + pow(x0, 3);
	zcfs[0] = d * _coeffs[0] + (b1 + b2) / 2.0;
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = _coeffs[i] * d;
	}
	rz = ISLIB_ABS(d)*_r + 0.5 * ISLIB_ABS(b1 - b2);
	return AAFR{ zidx, zcfs, rz };
}

//%============================================================================
//% End of module: revised affine forms
//%============================================================================