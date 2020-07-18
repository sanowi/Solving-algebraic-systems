//%
//%##########################################################################
//%
//%     Copyright (C) 2011 - Iwona Skalna
//%     AGH University of Science and Technology
//%     Department of Applied Computer Science
//%
//%     Class: Afine Forms 
//%            Based on papers of Comba & Stofli
//%
//%##########################################################################
//%

#include <cfloat>
#include <omp.h>
#include "../utils/stdafx.h"
#include "../utils/qsortm.h"
#include "multaux.h"
#include "aaf.h"

unsigned int AAF::_last = 0; //% initialization of the last index

int MULTS = MULTS_BEST; //% choice of multiplications formula

#define REDUCE //% if defined zeros are reduced in affine forms and in revised affine forms

const real zveceps = 1.0e-3; //% defined for Hladik's multplication

AAF::AAF(const interval& x)
//%----------------------------------------------------------------------------
//% constructor: affine form from interval
//%----------------------------------------------------------------------------
{
	if (x.rad() == 0) { //% x is a real number
		_indexes = cvector{ 1 };
		_coeffs = dvector{ 1 };
		_indexes[0] = 0;
		_coeffs[0] = x.mid();
	}
	else { //% x is a proper interval
		_indexes = cvector{ 2 };
		_coeffs = dvector{ 2 };
		_indexes[0] = 0;
		_indexes[1] = inclast();
		_coeffs[0] = x.mid();
		_coeffs[1] = x.rad();
	}
}

AAF::AAF(const AAF& x)
//%----------------------------------------------------------------------------
//% copying constructor
//%----------------------------------------------------------------------------
{
	_indexes = x._indexes;
	_coeffs = x._coeffs;
}

AAF&
//%----------------------------------------------------------------------------
//% copying operator
//%----------------------------------------------------------------------------
AAF::operator = (const AAF& x)
{
	if (this != &x) {
		_indexes = x._indexes;
		_coeffs = x._coeffs;
	}
	return (*this);
}

real
AAF::rad()
//%----------------------------------------------------------------------------
//% radius of affine form
//%----------------------------------------------------------------------------
{
	int	ltemp = _indexes.size();
	real r = 0.0;

	RoundUp();
	for (int i = 1; i < ltemp; ++i) {
		if (_coeffs[i] != 0.0) {
			r += ISLIB_ABS(_coeffs[i]);
		}
	}
	RoundNear();
	return r;
}

interval
//%----------------------------------------------------------------------------
//% Conversion: affine form => interval
//%----------------------------------------------------------------------------
AAF::reduce() const
{
	int	ltemp = _indexes.size();
	real r = 0.0, inf, sup;

	RoundUp(); //% set upward rounding mode
	for (int i = 1; i < ltemp; ++i) {
		if (_coeffs[i] != 0.0) {
			r += ISLIB_ABS(_coeffs[i]);
		}
	}
	inf = -(-_coeffs[0] + r);
	sup = _coeffs[0] + r;
	RoundNear(); //% set to nearest rounding mode 
	return interval{ inf, sup };
}

void
//%----------------------------------------------------------------------------
//% Eliminates zero coefficients
//%----------------------------------------------------------------------------
AAF::reduceZeroes()
{
	int	n = _indexes.size(), s = 0;
	cvector	_newIndexes;
	dvector	_newCoeffs;

	for (int i = 1; i < n; ++i) {
		if (_coeffs[i] == 0.0) ++s;
	}
	if (s == 0) return; //% there is no zero coefficients

	//% the zero coefficients should be removed from an affine form
	_newIndexes = cvector{ n - s };
	_newCoeffs = dvector{ n - s };
	_newIndexes.fill_in(0);
	_newCoeffs.fill_in(0.0);

	_newIndexes[0] = _indexes[0];
	_newCoeffs[0] = _coeffs[0];
	for (int i = 1, j = 1; i < n; ++i) {
		if (_coeffs[i] != 0.0) {
			_newIndexes[j] = _indexes[i];
			_newCoeffs[j] = _coeffs[i];
			++j;
		}
	}
	_indexes = _newIndexes;
	_coeffs = _newCoeffs;
}

void
AAF::setZeroAt(int i)
{
	_coeffs[i] = 0.0;
}

cvector
//%----------------------------------------------------------------------------
//% union of two sets of indexes
//%----------------------------------------------------------------------------
AAF::idx_union(cvector& a, cvector& b)
{
	int	ia = 0, ib = 0, ic = 0;
	cvector	c{ a.size() + b.size() };

	for (;;) {
		if (a[ia] == b[ib]) {
			c[ic++] = a[ia]; ia++; ib++;
		}
		else if (a[ia] < b[ib]) {
			c[ic++] = a[ia]; ia++;
		}
		else {
			c[ic++] = b[ib]; ib++;
		}
		if (ia >= a.size()) {
			while (ib < b.size()) {
				c[ic++] = b[ib]; ib++;
			}
			break;
		}
		if (ib >= b.size()) {
			while (ia < a.size()) {
				c[ic++] = a[ia]; ia++;
			}
			break;
		}
	}
	cvector wc{ ic };
	for (ia = 0; ia < ic; ia++) {
		wc[ia] = c[ia];
	}
	return wc;
}

dvector
//%----------------------------------------------------------------------------
//% Union of two sets of coefficients
//%----------------------------------------------------------------------------
AAF::cfs_union(AAF& a, AAF& b)
{
	int	ia = 0, ib = 0, ic = 0;
	dvector	c{ a._indexes.size() + b._indexes.size() };

	for (;;) {
		if (a._indexes[ia] == b._indexes[ib]) {
			c[ic++] = a._coeffs[ia++] + b._coeffs[ib++];
		}
		else if (a._indexes[ia] < b._indexes[ib]) {
			c[ic++] = a._coeffs[ia++];
		}
		else {
			c[ic++] = b._coeffs[ib++];
		}
		if (ia >= a._indexes.size()) {
			while (ib < b._indexes.size()) {
				c[ic++] = b._coeffs[ib++];
			}
			break;
		}
		if (ib >= b._indexes.size()) {
			while (ia < a._indexes.size()) {
				c[ic++] = a._coeffs[ia++];
			}
			break;
		}
	}
	dvector w{ ic };
	for (ia = 0; ia < ic; ia++) {
		w[ia] = c[ia];
	}
	return w;
}

AAF
//%----------------------------------------------------------------------------
//% Addition: affine form + affine form
//%----------------------------------------------------------------------------
operator + (const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real d = 0.0;
	cvector xidx{ x.idx() }, yidx{ y.idx() }, zidx;
	dvector	xcfs{ x.cfs() }, ycfs{ y.cfs() };

	zidx = AAF::idx_union(xidx, yidx); //% union of indexes
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }; //% nzcfs - near mode
	dvector uzcfs{ ltemp }; //% uzcfs - round-up mode
	dvector dzcfs{ ltemp }; //% dzcfs - round-down mode

	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = ycfs[iy++];  //% y.coeffs[ib]+0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++];  //% x.coeffs[ia]+0
			continue;
		}
		zcfs[iz] = xcfs[ix++] + ycfs[iy++];
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
		dzcfs[iz] = -(-xcfs[ix] - ycfs[iy]);
		uzcfs[iz] = xcfs[ix++] + ycfs[iy++];
		d += ISLIB_MAX(uzcfs[iz] - zcfs[iz], zcfs[iz] - dzcfs[iz]);
	}
	RoundNear();
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index used so far
		zcfs = zcfs.resize(1, d);            //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Addition: affine form + real number
//%----------------------------------------------------------------------------
operator + (const AAF& x, const real y)
{
	real zx, dz, uz, d = 0.0;
	cvector	zidx = x._indexes;
	dvector	zcfs = x._coeffs;

	zx = zcfs[0] + y;
	RoundUp();
	dz = -(-zcfs[0] - y);
	uz = zcfs[0] + y;
	d += ISLIB_MAX(uz - zx, zx - dz);
	RoundNear();
	zcfs[0] = zx;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index used so far
		zcfs = zcfs.resize(1, d);            //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Addition: real number + affine form
//%----------------------------------------------------------------------------
operator + (const real y, const AAF& x)
{
	real zx, dz, uz, d = 0.0;
	cvector	zidx = x._indexes;
	dvector	zcfs = x._coeffs;

	zx = y + zcfs[0];
	RoundUp();
	uz = y + zcfs[0];
	dz = -(-y - zcfs[0]);
	d += ISLIB_MAX(uz - zx, zx - dz);
	RoundNear();
	zcfs[0] = zx;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index used so far
		zcfs = zcfs.resize(1, d);            //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

void
//%----------------------------------------------------------------------------
//% Increment by affine form
//%----------------------------------------------------------------------------
operator += (AAF& x, const AAF& y)
{
	x = x + y;
}

void
//%----------------------------------------------------------------------------
//% Increment by real number
//%----------------------------------------------------------------------------
operator += (AAF& x, const real y)
{
	x = x + y;
}

AAF
//%----------------------------------------------------------------------------
//% Negation
//%----------------------------------------------------------------------------
operator - (const AAF& x)
{
	int	n = x._indexes.size();
	dvector	ncoeffs{ n };

	for (int i = 0; i < n; ++i) {
		ncoeffs[i] = -(x._coeffs[i]);
	}
	AAF z{ x._indexes, ncoeffs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Subtraction: affine form - affine form
//%----------------------------------------------------------------------------
operator - (const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x._indexes.size(), ny = y._indexes.size();
	real d = 0.0; //% error
	cvector xidx = x._indexes, yidx = y._indexes, zidx;
	dvector	xcfs = x._coeffs, ycfs = y._coeffs;

	zidx = AAF::idx_union(xidx, yidx); //% sum of indices
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }, uzcfs{ ltemp }, dzcfs{ ltemp };

	for (ix = 0, iy = 0, iz = 0; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = -ycfs[iy++];  //% 0 - y.coeffs[ib]
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++];  //% x.coeffs[ia] - 0
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
	RoundNear();
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last);  //% last index used so far
		zcfs = zcfs.resize(1, d);             //% rounding error       
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Subtraction: affine form - real number
//%----------------------------------------------------------------------------
operator - (const AAF& x, const real y)
{
	real zx, dz, uz, d = 0.0;
	cvector	zidx = x._indexes;
	dvector	zcfs = x._coeffs;

	zx = zcfs[0] - y;
	RoundUp();
	dz = -(-zcfs[0] + y);
	uz = zcfs[0] - y;
	d += ISLIB_MAX(uz - zx, zx - dz);
	RoundNear();
	zcfs[0] = zx;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index used so far
		zcfs = zcfs.resize(1, d);            //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Subtraction: real number - affine form
//%----------------------------------------------------------------------------
operator - (const real y, const AAF& x)
{
	real zx, dz, uz, d = 0.0;
	cvector	zidx = x._indexes;
	dvector	zcfs = x._coeffs;

	zx = y - zcfs[0];
	RoundUp();
	dz = -(-y + zcfs[0]);
	uz = y - zcfs[0];
	d += ISLIB_MAX(uz - zx, zx - dz);
	RoundNear();
	zcfs[0] = zx;
	for (int i = 1; i < zcfs.size(); ++i) {
		zcfs[i] = -zcfs[i];
	}
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index used so far
		zcfs = zcfs.resize(1, d);            //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Verified trivial multiplication of affine forms.
//%----------------------------------------------------------------------------
AAF::quickMult(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real d, sx, sy, sx2 = 0.0, sy2 = 0.0;
	cvector	xidx = x.idx(), yidx = y.idx(), zidx;
	dvector	xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp }, xcfs2{ ltemp }, ycfs2{ ltemp };

	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	//% create new arrays so that the all the indices from 1 to n are present
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) { //% x index does not exist
			ycfs2[iz] = ycfs[iy++];
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) { // y% index does not exist
			xcfs2[iz] = xcfs[ix++];
			continue;
		}
		xcfs2[iz] = xcfs[ix++];
		ycfs2[iz] = ycfs[iy++];
	}
	//% RoundNear mode
	zcfs[0] = xcfs[0] * ycfs[0];
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	sx = sy = d = 0.0;
	RoundUp(); //% ROundUp mode
	dzcfs[0] = -(-xcfs[0] * ycfs[0]); //% minus added to obtain round down mode
	uzcfs[0] = xcfs[0] * ycfs[0]; //% round up
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]); //% partial round-off error
	for (int i = 1; i < ltemp; ++i) {
		dzcfs[i] = -(-xcfs[0] * ycfs2[i] + (-ycfs[0] * xcfs2[i])); //% round down
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i]; //% round up
		d += ISLIB_MAX(uzcfs[i] - zcfs[i], zcfs[i] - dzcfs[i]); //% partial round-off error
		sx += ISLIB_ABS(xcfs2[i]);
		sy += ISLIB_ABS(ycfs2[i]);
	}
	d += sx * sy;
	RoundNear();
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last free index
		zcfs = zcfs.resize(1, d);            //% final round off error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Verified multiplication of affine forms based on Miyaima's et al result 
//% from paper "On the best multiplication of affine forms"
//%----------------------------------------------------------------------------
AAF::bestMult(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real ez, dd, dmin, dmax;
	interval d, rx, ry;
	cvector	xidx = x.idx(), yidx = y.idx(), zidx;
	dvector	xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector	xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp };

	//% the case when one of the factors is in fact a real number
	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) return (rx.mid() * ry.mid());
	if (rx.rad() == 0.0) return rx.mid() * y;
	if (ry.rad() == 0.0) return x * ry.mid();

	//% the case when both affine forms are non-trivial
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	//% create new arrays so that all the indices from 1 to n are present
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
	dMinMax(0.0, 0.0, xcfs2, ycfs2, dmin, dmax); //% this gives minimal range
	d = interval{ dmin, dmax };
	RoundNear(); //% is this OK?
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	dd = 0.0;
	RoundUp();
	dzcfs[0] = -(-xcfs[0] * ycfs[0] - (d.inf() - (-d.sup() + d.inf()) / 2.0));
	uzcfs[0] = xcfs[0] * ycfs[0] + d.inf() + (d.sup() - d.inf()) / 2.0;
	dd += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (int i = 1; i < ltemp; ++i) {
		dzcfs[i] = -(-xcfs[0] * ycfs2[i] - ycfs[0] * xcfs2[i]);
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
		dd += ISLIB_MAX(uzcfs[i] - zcfs[i], zcfs[i] - dzcfs[i]);
	}
	ez = dd + d.rad();
	RoundNear();
	if (ez != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% ++AAF::_last is the last unused index
		zcfs = zcfs.resize(1, ez);           //% noise variable with rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Verified standard multiplication of affine forms (Kolev's version)
//%----------------------------------------------------------------------------
AAF::stdMult(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real d = 0.0, s = 0.0, s2 = 0.0, s3 = 0.0, sup = 0.0, sdown = 0.0, sabs = 0.0, es, ds;
	cvector xidx = x.idx(), yidx = y.idx(), zidx;
	dvector xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx);
	ltemp = zidx.size();
	dvector	zcfs{ ltemp }, uzcfs{ ltemp }, dzcfs{ ltemp };

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum(x_i * y_i)
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * s;
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		sdown += xcfs[ix] * ycfs[iy];
		sabs += ISLIB_ABS(xcfs[ix] * ycfs[iy]);
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sdown;

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	es = 0.5 * sabs; //% ISLIB_ABS(sdown);
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]);
			iy++;  continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]);
			ix++;  continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); iy++;
		s2 += ISLIB_ABS(xcfs[ix]); ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sup;
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
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
	ds = s2 * s3 + d;
	RoundNear(); //% ???
	ds -= es;
	if (ds != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index use so far
		zcfs = zcfs.resize(1, ds);           //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Multiplication of affine forms based on tangent plane
//%----------------------------------------------------------------------------
AAF::AASEEMult(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real d = 0.0, s = 0.0, s2 = 0.0, s3 = 0.0, sup = 0.0, sdown = 0.0, sabs = 0.0, ds;
	cvector xidx = x.idx(), yidx = y.idx(), zidx;
	dvector xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx);
	ltemp = zidx.size();

	dvector	zcfs{ ltemp }, uzcfs{ ltemp }, dzcfs{ ltemp };

	//% Change: first check ltemp i and return the form immediately
	if (ltemp == 2) {
		return stdMult(x, y);
	}
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum(x_i * y_i)
	}
	zcfs[0] = xcfs[0] * ycfs[0] + s;

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		sdown += xcfs[ix] * ycfs[iy];
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + sdown;

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;  continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;  continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		ix++; iy++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + sup;

	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
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
	real s4 = 0.0;
	//% quadratic time complexity
	for (int i = 1; i < ltemp; ++i) {
		for (int j = 1; j < i; ++j) {
			if (i != j) {
				s4 += ISLIB_ABS(xcfs2[i] * ycfs2[j] + xcfs2[j] * ycfs2[i]);
			}
		}
	}
	ds = s4 + d;

	RoundNear(); //% ???
	if (ds != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index use so far
		zcfs = zcfs.resize(1, ds);           //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Verified standard multiplication of affine forms (Kolev's version).
//% Improved???
//%----------------------------------------------------------------------------
AAF::stdMultImp(const AAF& x, const AAF& y)
{
	int ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real d = 0.0, s = 0.0, s2 = 0.0, s3 = 0.0, sup = 0.0, sdown = 0.0, sabs = 0.0, es, ds;
	cvector	 xidx = x.idx(), yidx = y.idx(), zidx;
	dvector xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx);
	ltemp = zidx.size();
	dvector zcfs{ ltemp }, uzcfs{ ltemp }, dzcfs{ ltemp };

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++;
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++;
			continue;
		}
		s += xcfs[ix++] * ycfs[iy++]; //% sum(x_i * y_i)
	}
	zcfs[0] = xcfs[0] * ycfs[0] + 0.5 * s;
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			zcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			zcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		zcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	RoundDown();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			iy++; continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			ix++; continue;
		}
		sdown += xcfs[ix] * ycfs[iy];
		ix++;
		iy++;
	}
	dzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sdown;

	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			dzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
			continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			dzcfs[iz] = xcfs[ix++] * ycfs[0];  //% x.coeffs[ ia ] * 0
			continue;
		}
		dzcfs[iz] = xcfs[ix++] * ycfs[0] + xcfs[0] * ycfs[iy++];
	}
	es = 0.5 * ISLIB_ABS(sdown);
	RoundUp();
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			s3 += ISLIB_ABS(ycfs[iy]);
			iy++;  continue;
		}
		if (iy == ny || yidx[iy] != zidx[iz]) {
			s2 += ISLIB_ABS(xcfs[ix]);
			ix++;  continue;
		}
		sup += xcfs[ix] * ycfs[iy];
		s3 += ISLIB_ABS(ycfs[iy]); iy++;
		s2 += ISLIB_ABS(xcfs[ix]); ix++;
	}
	uzcfs[0] = xcfs[0] * ycfs[0] + 0.5 * sup;
	d += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != zidx[iz]) {
			uzcfs[iz] = xcfs[0] * ycfs[iy++];  //% y.coeffs[ ib ] * 0 <=> x.coeff[ ia ] = 0
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
	ds = s2 * s3 + d;
	RoundNear(); //% ???
	ds -= es;
	if (ds != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% last index use so far
		zcfs = zcfs.resize(1, ds);           //% rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Hladik's multiplication
//%----------------------------------------------------------------------------
AAF::hlMult(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real ez, dd;
	interval d, rx, ry;
	cvector	xidx = x.idx(), yidx = y.idx(), zidx;
	dvector	xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector	xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp };

	//% the case when one of the factors is in fact a real number
	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) return (rx.mid() * ry.mid());
	if (rx.rad() == 0.0) return rx.mid() * y;
	if (ry.rad() == 0.0) return x * ry.mid();

	//% the case when both affine forms are non-trivial
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	//% create new arrays so that all the indices from 1 to n are present
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
	//% the vectors xcfs2 and ycfs2 have the same length
	//% n/2[x^Ty-||x||*||y||,x^y+||x||*||y||
	dvector vx(ltemp - 1), vy(ltemp - 1);
	for (int i = 1; i < ltemp; ++i) {
		vx[i - 1] = xcfs2[i];
		vy[i - 1] = ycfs2[i];
	}
	real xty = 0.0;
	real xnorm = norm(vx);
	real ynorm = norm(vy);

	for (int i = 1; i < ltemp; ++i) {
		xty += xcfs2[i] * ycfs2[i];
	}
	real d0 = xnorm * ynorm;
	d = (0.5 * (ltemp - 1)) * interval{ xty - d0, xty + d0 };
	RoundNear(); //% ???
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	dd = 0.0;
	RoundUp();
	dzcfs[0] = -(-xcfs[0] * ycfs[0] - (d.inf() - (-d.sup() + d.inf()) / 2.0));
	uzcfs[0] = xcfs[0] * ycfs[0] + d.inf() + (d.sup() - d.inf()) / 2.0;
	dd += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (int i = 1; i < ltemp; ++i) {
		dzcfs[i] = -(-xcfs[0] * ycfs2[i] - ycfs[0] * xcfs2[i]);
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
		dd += ISLIB_MAX(uzcfs[i] - zcfs[i], zcfs[i] - dzcfs[i]);
	}
	ez = dd + d.rad();
	RoundNear();
	if (ez != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% ++AAF::_last is the last unused index
		zcfs = zcfs.resize(1, ez);           //% noise variable with rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

void
dhlMS(int ltemp2, const dvector& vx, const dvector& vy, const dvector& zvec, interval& d)
{
	real xtzzy = 0.0, d0;

	dvector zx{ ltemp2 }, zy{ ltemp2 };
	for (int i = 0; i < ltemp2; ++i) {
		zx[i] = vx[i] / zvec[i];
		zy[i] = vy[i] / zvec[i];
	}

	real ztz = 0.0, xzzy = 0.0;
	for (int i = 0; i < ltemp2; ++i) {
		ztz += zvec[i] * zvec[i]; //% z^T*z
		xzzy += zx[i] * zy[i];
	}
	d0 = norm(zx)*norm(zy); //% ||Z^{-1}x||*||Z^{-1}y||
	d = 0.5 * ztz * interval{ xzzy - d0, xzzy + d0 };
}

AAF
//%----------------------------------------------------------------------------
//% Hladik's multiplication
//%----------------------------------------------------------------------------
AAF::hlMultScaled(const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x.idx().size(), ny = y.idx().size();
	real ez, dd;
	interval d, d2, rx, ry;
	cvector	xidx = x.idx(), yidx = y.idx(), zidx;
	dvector	xcfs = x.cfs(), ycfs = y.cfs();

	zidx = AAF::idx_union(xidx, yidx); //% sum of affine forms indices
	ltemp = zidx.size();

	dvector	xcfs2{ ltemp }, ycfs2{ ltemp }, zcfs{ ltemp }, dzcfs{ ltemp }, uzcfs{ ltemp };

	//% the case when one of the factors is in fact a real number
	rx = x.reduce();
	ry = y.reduce();
	if (rx.rad() == 0.0 && ry.rad() == 0.0) return (rx.mid() * ry.mid());
	if (rx.rad() == 0.0) return rx.mid() * y;
	if (ry.rad() == 0.0) return x * ry.mid();

	//% the case when both affine forms are non-trivial
	xcfs2.fill_in(0.0);
	ycfs2.fill_in(0.0);
	//% create new arrays so that all the indices from 1 to n are present
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
	//% the vectors xcfs2 and ycfs2 have the same length

	//% n/2[x^Ty-||x||*||y||,x^y+||x||*||y||
	int ltemp2 = ltemp - 1;
	dvector vx(ltemp2), vy(ltemp2);
	for (int i = 1; i < ltemp; ++i) {
		vx[i - 1] = xcfs2[i];
		vy[i - 1] = ycfs2[i];
	}

	dvector zmin{ ltemp2 }, zmax{ ltemp2 };
	//real xnorm, ynorm;
	// computing eigenvectors:
	// "for" min : zmin = ||y||x-||x||y
	// "for" max : zmax = ||y||x+||x||y
	/*xnorm = norm(vx); ynorm = norm(vy);
	for (int i = 0; i < ltemp2; ++i) {
		zmin[i] = ynorm * vx[i] - xnorm * vy[i];
		if (zmin[i] < 0.0001) zmin[i] = zveceps;
		zmax[i] = ynorm * vx[i] + xnorm * vy[i];
		if (zmax[i] < 0.0001) zmax[i] = zveceps;
	}*/
	/*std::cout << "zmin = " << zmin << std::endl;
	std::cout << "zmax = " << zmax << std::endl;
	char c = _getch();*/
	zmin.fill_in(10); zmax.fill_in(10);

	dhlMS(ltemp2, vx, vy, zmin, d2);
	dhlMS(ltemp2, vx, vy, zmax, d);
	d = interval{ ISLIB_MAX(d.inf(), d2.inf()), ISLIB_MIN(d.sup(), d2.sup()) };

	RoundNear(); //% for now hereu, but perhaps should be placed somewhere else
	zcfs[0] = xcfs[0] * ycfs[0] + d.mid();
	for (int i = 1; i < ltemp; ++i) {
		zcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
	}
	dd = 0.0;
	RoundUp();
	dzcfs[0] = -(-xcfs[0] * ycfs[0] - (d.inf() - (-d.sup() + d.inf()) / 2.0));
	uzcfs[0] = xcfs[0] * ycfs[0] + d.inf() + (d.sup() - d.inf()) / 2.0; //% d.mid() = d.inf() + (d.sup() - d.inf())/2
	dd += ISLIB_MAX(uzcfs[0] - zcfs[0], zcfs[0] - dzcfs[0]);
	for (int i = 1; i < ltemp; ++i) {
		dzcfs[i] = -(-xcfs[0] * ycfs2[i] - ycfs[0] * xcfs2[i]);
		uzcfs[i] = xcfs[0] * ycfs2[i] + ycfs[0] * xcfs2[i];
		dd += ISLIB_MAX(uzcfs[i] - zcfs[i], zcfs[i] - dzcfs[i]);
	}
	ez = dd + d.rad();
	RoundNear();
	if (ez != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last); //% ++AAF::_last is the last unused index
		zcfs = zcfs.resize(1, ez);           //% noise variable with rounding error
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Multiplication of affine forms
//% could be: quick, best, standard, minrange
//%----------------------------------------------------------------------------
operator * (const AAF& x, const AAF& y)
{
	switch (MULTS) {
	case MULTS_QUICK: // QUICK MULT
		return AAF::quickMult(x, y);
	case MULTS_BEST: // BEST MULT
		return AAF::bestMult(x, y);
	case MULTS_STD: // STD MULT
		return AAF::stdMult(x, y);
	case MULTS_MH: // ShlTD MULT
		return AAF::hlMult(x, y);
	case MULTS_MHS: // ShlTD MULT
		return AAF::hlMultScaled(x, y);
	case MULTS_AASEE: // tangent plane, MULTS_AASEE
		return AAF::AASEEMult(x, y);
	default:
		return AAF::bestMult(x, y);
	}
}

AAF
//%----------------------------------------------------------------------------
//% Multiplication of an affine form and a real number (on the left)
//% with rounging errors
//%----------------------------------------------------------------------------
operator * (const real x, const AAF& y)
{
	int	ltemp = y._indexes.size();
	real d = 0.0;
	cvector zidx = y._indexes;
	dvector	zcfs = y._coeffs, zx{ ltemp }, dz{ ltemp }, uz{ ltemp };
	
	for (int i = 0; i < ltemp; ++i) {
		zx[i] = zcfs[i] * x;
	}
	RoundUp();
	for (int i = 0; i < ltemp; ++i) {
		dz[i] = -zcfs[i] * x;
		uz[i] = zcfs[i] * x;
		d += ISLIB_MAX(uz[i] - zx[i], zx[i] + dz[i]);
		zcfs[i] = zx[i];
	}
	RoundNear();
	if (d != 0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Multiplication of an affine form and a real number (on the right)
//% with rounging errors
//%----------------------------------------------------------------------------
operator * (const AAF& y, const real x)
{
	int ltemp = y._indexes.size();
	real d = 0.0;
	cvector zidx = y._indexes;
	dvector zcfs = y._coeffs, zx{ ltemp }, dz{ ltemp }, uz{ ltemp };
	
	for (int i = 0; i < ltemp; ++i) {
		zx[i] = zcfs[i] * x;
	}
	RoundUp();
	for (int i = 0; i < ltemp; ++i) {
		dz[i] = -zcfs[i] * x;
		uz[i] = zcfs[i] * x;
		d += ISLIB_MAX(uz[i] - zx[i], zx[i] + dz[i]);
		zcfs[i] = zx[i];
	}
	RoundNear();
	if (d != 0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Chebyshev approximation of f( x ) = 1 / x reciprocal of affine form 
//% (Kolev's, Miyajima's version) approximation: ax + b + d*eps_new
//% calculated as an average of tangent and secant
//% with rounging errors
//%----------------------------------------------------------------------------
AAF::reciprocal() const
{
	int	i, ltemp = _indexes.size();
	interval c = reduce(); //% range of affine form
	real a, x, x0, b1, b2, b, rw;
	cvector	ridx = _indexes; //% indexes of affine form
	dvector	rcfs{ ltemp };
	
	a = -1.0 / (c.inf() * c.sup()); //% 1/(\u x*\o x), gradient
	x = -sqrt(-1.0 / a);
	x0 = (c > 0.0) ? -x : x; //% x0 - tangent point
	b1 = 1.0 / c.inf() - a * c.inf(); //% b1=1/(\u x)+1/(\o x)
	b2 = 1.0 / x0 - a * x0; //% b2=1/sqrt(\u x*\o x)+1/(\u x*\o x)*sqrt(\u x*\o x)
	b = b2 + 0.5 * (b1 - b2); //% 0.5*(b2+b1);
	rw = b - b2; //% b-b2;
	rcfs[0] = a * (_coeffs[0]) + b;
	for (i = 1; i < ltemp; ++i) {
		rcfs[i] = a * (_coeffs[i]);
	}
	if (rw != 0.0) {
		ridx = ridx.resize(1, ++AAF::_last);
		rcfs = rcfs.resize(1, rw);
	}
	AAF z{ ridx, rcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Min-range approximation of f( x ) = 1 / x
//% reciprocal of affine form (Kolev's, Miyajima's version)
//% approximation: ax + b + d*eps_new (Chebyshev approximation)
//% calculated as average of secant and tangent
//%----------------------------------------------------------------------------
AAF::reciprocal2() const
{
	int	i, ltemp = _indexes.size();
	interval c = reduce();
	real xd, xu, a, b1, b2, rw;
	cvector	ridx = _indexes;
	dvector	rcfs{ ltemp };

	//RoundUp();
	xd = c.inf(); //% \u(x)
	xu = c.sup(); //% \o(x)
	b1 = xu * xu;
	b2 = 2.0 * xd * b1; //% 2*\u(y)*\o(y)^2
	a = -1.0 / b1;
	rcfs[0] = a * (_coeffs[0]) + (xd + xu) * (xd + xu) / b2;
	rw = (xd - xu) * (xd - xu) / b2;
	for (i = 1; i < ltemp; ++i) {
		rcfs[i] = a * (_coeffs[i]);
	}
	RoundNear();
	if (rw != 0.0) {
		ridx = ridx.resize(1, ++AAF::_last);
		rcfs = rcfs.resize(1, rw);
	}
	AAF z{ ridx, rcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Division of affine forms
//% Kolev, Hasan El Owny (patrz El-Owny) version
//% round off errors not included yet
//%----------------------------------------------------------------------------
operator / (const AAF& x, const AAF& y)
{
	int	ix, iy, iz, ltemp, nx = x._indexes.size(), ny = y._indexes.size();
	real c;
	cvector	xidx = x._indexes, yidx = y._indexes, zidx, pidx;
	dvector xcfs = x._coeffs, ycfs = y._coeffs, pcfs;
	AAF p{ 0 }, ry{ 0 }, pry{ 0 };

	pidx = AAF::idx_union(xidx, yidx); // sum of indexes
	ltemp = pidx.size();
	pcfs = dvector{ ltemp };
	c = (x._coeffs[0]) / (y._coeffs[0]); // y._coeffs[0] cannot be zero, since [y] does not contain zero
	pcfs[0] = 0.0;
	for (ix = 1, iy = 1, iz = 1; iz < ltemp; ++iz) {
		if (ix == nx || xidx[ix] != pidx[iz]) {
			pcfs[iz] = -(c * ycfs[iy++]);
			continue;
		}
		if (iy == ny || yidx[iy] != pidx[iz]) {
			pcfs[iz] = xcfs[ix++];
			continue;
		}
		pcfs[iz] = xcfs[ix++] - (c * ycfs[iy++]);
	}
	p = AAF{ pidx, pcfs };
	ry = y.reciprocal(); // or can be reciprocal (Chebyshev), reciprocal2 (min-range)
	pry = p * ry;
	pry._coeffs[0] += c;
#ifdef REDUCE
	pry.reduceZeroes();
#endif
	return pry;
}

//%============================================================================
//% FUNCTIONS on affine forms
//%----------------------------------------------------------------------------
//%
//% * sqrtf - square root (Chebyshev approximation)
//% * sqrttf - square root (Min-range approximation)
//% * sqrf - square function (Chebyshev approximation)
//% * sqrrf - square function (Min-range approximation)
//% * exp - exponent (Chebyshev approximation)
//%
//%============================================================================

AAF
//%----------------------------------------------------------------------------
//% Chebyshev approximation of square root of an affine form
//% algorithm (a) from AA2_Wu_Hroud paper
//%----------------------------------------------------------------------------
AAF::sqrtf()
{
	int	i, ltemp = _indexes.size();
	real a, b, alpha, beta, delta, db, dmin, dmax, fad, fau, fbd, fbu, aad, aau, abd, abu, gau, dfau;
	interval xi = reduce();
	cvector zidx = _indexes;
	dvector zcfs{ ltemp };

	a = xi.inf();
	b = xi.sup();
	RoundUp();
	fau = sqrt(a);
	fbu = sqrt(b);
	alpha = (sqrt(b) - sqrt(a)) / (b - a);
	aau = alpha * a;
	abu = alpha * b;
	gau = 1.0 / (4.0 * alpha);
	dfau = 1.0 / (2.0 * sqrt(a));
	RoundDown();
	fad = sqrt(a);
	fbd = sqrt(b);
	aad = alpha * a;
	abd = alpha * b;
	RoundNear();
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
	delta = interval{ dmin, dmax }.rad();
	for (i = 0; i < ltemp; ++i) {
		zcfs[i] = alpha * (_coeffs[i]);
	}
	zcfs[0] += beta;
	zidx = zidx.resize(1, ++AAF::_last);
	zcfs = zcfs.resize(1, delta);
	AAF z{ zidx, zcfs };
#ifdef REDUCE
	z.reduceZeroes();
#endif
	return z;
}

AAF
//%----------------------------------------------------------------------------
//% Min-range approximation of square root of an affine form
//%----------------------------------------------------------------------------
AAF::sqrttf()
{
	int i, ltemp = _indexes.size();
	real a, b, alpha, b1, b2, d;
	interval xi = reduce();
	cvector zidx = _indexes;
	dvector zcfs{ ltemp };

	a = xi.inf();
	b = xi.sup();
	//RoundUp(); //% ???
	alpha = 1.0 / (2.0 * sqrt(b));
	b1 = 0.5 * sqrt(b);
	b2 = sqrt(a) - (a * alpha);
	//RoundDown();  //% ???
	//RoundNear();  //% ???
	for (i = 0; i < ltemp; ++i) {
		zcfs[i] = alpha * (_coeffs[i]);
	}
	zcfs[0] += (b1 + b2) / 2.0;
	d = (b1 - b2) / 2.0;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	return AAF{ zidx, zcfs };
}

AAF
//%----------------------------------------------------------------------------
//% Chebyshev (min-max) approximation of square of an affine form
//%----------------------------------------------------------------------------
AAF::sqrf()
{
	int ltemp = _indexes.size();
	real cinf, csup, b1, b2, d, a;
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
		zcfs[i] = (_coeffs[i]) * a;
	}
	zcfs[0] += (b1 + b2) / 2.0;
	d = ISLIB_ABS(b1 - b2) / 2.0;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	return AAF{ zidx, zcfs };
}

AAF
//%----------------------------------------------------------------------------
//% Min-range approximation square of affine form
//%----------------------------------------------------------------------------
AAF::sqrrf()
{
	int	ltemp = _indexes.size();
	interval c = reduce();
	real cinf, csup, b1, b2, d, a;
	cvector	zidx = _indexes;
	dvector	zcfs{ ltemp };

	cinf = c.inf();
	csup = c.sup();
	a = 2.0 * cinf;
	b1 = -(cinf * cinf);
	b2 = csup * csup - a * csup;
	for (int i = 0; i < ltemp; ++i) {
		zcfs[i] = (_coeffs[i]) * a;
	}
	zcfs[0] += (b1 + b2) / 2.0;
	d = ISLIB_ABS(b1 - b2) / 2.0;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	return AAF{ zidx, zcfs };
}

AAF
//%----------------------------------------------------------------------------
//% Chebyshev approximation of exponent of affine form
//%----------------------------------------------------------------------------
AAF::expf()
{
	int	ltemp = _indexes.size();
	real cinf, csup, b1, b2, d, a, yinf, ysup;
	interval c = reduce();
	cvector	zidx = _indexes;
	dvector	zcfs{ ltemp };

	cinf = c.inf();
	csup = c.sup();
	yinf = exp(cinf);
	ysup = exp(csup);
	a = (ysup - yinf) / (csup - cinf);
	if (a == 0.0) {
		return AAF(yinf);
	}
	b1 = (yinf * csup - ysup * cinf) / (csup - cinf);
	b2 = a * (1.0 - log(a));
	for (int i = 0; i < ltemp; ++i) {
		zcfs[i] = (_coeffs[i]) * a;
	}
	zcfs[0] += (b1 + b2) / 2.0;
	d = ISLIB_ABS(b1 - b2) / 2.0;
	if (d != 0.0) {
		zidx = zidx.resize(1, ++AAF::_last);
		zcfs = zcfs.resize(1, d);
	}
	return AAF{ zidx, zcfs };
}

std::ostream&
//%----------------------------------------------------------------------------
//% Print out an affine form
//%----------------------------------------------------------------------------
operator << (std::ostream& os, const AAF& a)
{
	dvector cfs = a._coeffs;
	cvector idx = a._indexes;

	os << cfs[0];
	for (int i = 1; i < cfs.size(); ++i) {
		if (cfs[i] >= 0) {
			os << " + " << cfs[i];
		}
		else if (cfs[i] < 0) {
			os << " - " << (-cfs[i]);
		}
		os << "e" << (idx[i]);
	}
	return os;
}

//%============================================================================
//% End of module: standard affine forms
//%============================================================================