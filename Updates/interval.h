//##########################################################################
//%
//% Copyright (C) 2011 - Iwona Skalna
//% AGH University of Science and Technology
//% Department of Applied Computer Science
//%
//% islib - library for proper intervals
//% template of _interval class
//%
//##########################################################################

#pragma once

//% MSVC++ 9.0   _MSC_VER = 1500
//% MSVC++ 8.0   _MSC_VER = 1400
//% MSVC++ 7.1   _MSC_VER = 1310
//% MSVC++ 7.0   _MSC_VER = 1300
//% MSVC++ 6.0   _MSC_VER = 1200
//% MSVC++ 5.0   _MSC_VER = 1100

#ifdef _MSC_VER				//% if Visual C++ compiler
#pragma fenv_access (on)	//% disables (ON) or enables (OFF) optimizations 
#endif						//% that could change flag tests and mode changes

#include <iostream>
#include <float.h>
#include <cassert>
#include "Round.h"

typedef double real;

extern real vnEta;
extern real vnEps;

real	nEps();
real	nEta();
real	predR(const real);
real	succR(const real);

namespace islib_interval
{
	//%//////////////////////////////////////////////////////////////
	//% BIG_ENDIAN and LITTLE_ENDIAN are checked 
	//% to reverse the order of ulongs
	//%//////////////////////////////////////////////////////////////
	const unsigned char endian_test[2] = { 1, 0 };
	const unsigned short ISLIB_ENDIANTEST = *(unsigned short*)endian_test;
	const unsigned short ISLIB_ONE = (unsigned short)1;
	//%==============================================================
	//% Special values
	//% Not a Number (NAN)
	//% Infinity (INF)
	//---------------------------------------------------------------
#if ( ISLIB_ENDIANTEST == ISLIB_ONE )
	//% LITTLE ENDIAN - e.g., Windows
	const unsigned long islib_inf[2] = { 0x00000000, 0x7ff00000 };
	const unsigned long islib_nan[2] = { 0xffffffff, 0x7fffffff };
#else
	//% BIG ENDIAN - Linux (Unix) like systems
	const unsigned long islib_inf[2]={0x7ff00000, 0x00000000};
	const unsigned long islib_nan[2]={0x7fffffff, 0xffffffff};
#endif
	const real ISLIB_INFINITY = *(real*)islib_inf;
	const real ISLIB_NAN = *(real*)islib_nan;
	//%==============================================================

	//%//////////////////////////////////////////////////////////////

	//%==============================================================
	//% MAX(x, y)
	//%--------------------------------------------------------------
	template <typename T>
	inline T const& ISLIB_MAX(T const& x, T const& y)
	{
		return ((x) < (y)) ? (y) : (x);
	}
	//%==============================================================
	//% MIN(x, y)
	//%--------------------------------------------------------------
	template <typename T>
	inline T const& ISLIB_MIN(T const& x, T const& y)
	{
		return ((x) > (y)) ? (y) : (x);
	}
	//%=============================================================
	//% ABS(x, y) absolute value
	//%-------------------------------------------------------------
	inline double ISLIB_ABS(const real x) // ABS(x)
	{
		return ((x) >= 0) ? (x) : (-x);
	}

	//%//////////////////////////////////////////////////////////////
	//%
	//% Interval Class 
	//% Template
	//%
	//%//////////////////////////////////////////////////////////////
	template <class T>
	class _interval
	{
	public:
		//% constructors
		_interval()	{ _inf = 0.0; _sup = 0.0; } //% zero interval
		inline _interval(const T& a);
		inline _interval(const T& a, const T& b);
		template<class T1, class T2> _interval(const T1 &a, const T2 &b);

		//% template<_interval<class T1>(const T1& a, const T& b) { _inf = a; _sup = b; }

		//% class methods
		T		inf() const { return _inf; }	//% left endpoint
		T		sup() const { return _sup; }	//% right endpoint
		bool	intersects(const _interval&, _interval&) const;
		void	succ(_interval&) const;
		void	pred(_interval&) const;
		inline	T rad() const;			//% radius
		inline	T width() const;		//% width
		inline	T mid() const;			//% midpoint
		inline	T mag()const;			//% maximal absolute value
		inline	T mig() const;			//% minimal absolute value
		inline	T dist(const _interval&) const; //% Hausdorff distance

		_interval reciprocal() const;	//% reciprocal of an interval
		_interval operator - () const;	//% negation operator
		_interval operator ~ () const;	//% dual operator

		inline bool is_empty() const;	//% checks if an interval is empty
		inline bool is_proper() const;	//% checks if an interval is proper

		bool set_contains(const T&) const;	//% inclusion
		bool set_contains(const _interval& T) const; //% set inclusion
		bool set_strictly_contains(const T&) const; //% strict inclusion
		bool set_strictly_contains(const _interval& T) const; //% set strict inclusion
		static const _interval emptyset() { return _interval(ISLIB_NAN); } //% empty set
		static const _interval universe() { return _interval(-ISLIB_INFINITY, ISLIB_INFINITY); } //% universe
		//template<class S> operator _interval<S>() { return _interval<S>((S) _inf, (S) _sup); } //% cast operator
	protected:
		T _inf;
		T _sup;
	};

	template<class T>
	inline _interval<T>::_interval(const T& a)
	{
#ifdef ISLIB_PRED
		_inf = predR(a); _sup = succR(a); 
#else
		_inf = a; _sup = a;
#endif
	}

	template<class T>
	inline _interval<T>::_interval(const T& a, const T& b)
	{
#ifdef ISLIB_PRED
		_inf = predR(a); _sup = succR(b);
#else
		_inf = a; _sup = b;
#endif
	}

	template<class T> template<class T1, class T2>
	_interval<T>::_interval(const T1 &a, const T2 &b)
	{
#ifdef ISLIB_PRED
		_inf = predR(a); _sup = succR(b);
#else
		_inf = a; _sup = b;
#endif
	}

	template<class T>
	bool _interval<T>::intersects(const _interval<T>& a, _interval<T>& x) const
	// intersection of [this] and [a]: [x] = [this] n [a]
	{
		real b = ISLIB_MAX(_inf, a._inf);
		real c = ISLIB_MIN(_sup, a._sup);
		if (b <= c) {
			x = _interval{ b, c };	// intersection of [this] and [x]
			return true;			// they intersect
		}
		x = emptyset();
		return false;				// they don't intersect, 
		// in other words the intersection is empty
	}

	template<class T>
	void _interval<T>::pred(_interval<T>& a) const
	// R = largest possible interval with R.inf > A.inf and R.sup < A.sup
	// Result: 1, if R is a true interval, 0 otherwise
	{
		real vnEta = ulp(1.0);
		RoundUp();
		real x = -(-_inf - vnEta);
		real y = _sup - vnEta;
		a = interval{ x, y };
		RoundNear();
		return (a._inf <= a._sup);
	}

	template<class T>
	void _interval<T>::succ(_interval<T>& a) const
	// R = smallest possible interval with R.inf < A.inf and R.sup > A.sup
	{
		real vnEta = ulp(1.0);
		RoundUp();
		real x = -(-_inf + vnEta);
		real y = _sup + vnEta;
		RoundNear();
		a = interval{ x, y };
	}

	template <class T>
	T _interval<T>::rad() const
	// radius of an interval
	{
		if (is_empty()) {
			return -1.0;
		}
		RoundUp();
		T r = ISLIB_ABS((_inf + (_sup - _inf) / (T)2) - _inf);
		// r = abs(_sup - _inf)/2; - the old version of midpoint
		RoundNear();
		return r;
	}

	template <class T>
	T _interval<T>::width() const
	// width of an interval
	{
		if (is_empty()) {
			return -1.0;
		}
		RoundUp();
		T r = _sup - _inf;
		RoundNear();
		return r;
	}

	template <class T>
	T _interval<T>::mid() const
	// midpoint of an interval suggested by Neumaier
	{
		if (is_empty()) {
			return ISLIB_NAN;
		}
		RoundUp();
		T m = _inf + (_sup - _inf) / (T)2;
		RoundNear();
		return m;
	}

	template <class T>
	T _interval<T>::mag() const
	// magnitude = max{ |x| : x in[x]}
	{
		if (is_empty()) {
			return ISLIB_NAN;
		}
		return ISLIB_MAX(ISLIB_ABS(_inf), ISLIB_ABS(_sup));
	}

	template <class T>
	T _interval<T>::mig() const
	// mignitude = min{ |x| : x in[x]} if 0 notin [x], 0 otherwise
	{
		if (is_empty()) {
			return ISLIB_NAN;
		}
		return (set_contains(0)) ? 0 : ISLIB_MIN(ISLIB_ABS(_inf), ISLIB_ABS(_sup));
	}

	template <class T>
	T _interval<T>::dist(const _interval<T>& x) const
	// Hausdorff distance
	{
		return ISLIB_MAX(ISLIB_ABS(_inf - x._inf), ISLIB_ABS(_sup - x._sup));
	}

	template <class T>
	bool _interval<T>::is_empty() const
	// checks if an interval is empty
	{
		return (_isnan(_inf) && _isnan(_sup) || _sup < _inf);
	}

	template <class T>
	bool _interval<T>::is_proper() const
		/* checks if an interval is proper */
	{
		if (is_empty())
		{
			return true;
		}
		return (_inf <= _sup);
	}

	template <class T>
	void hull(_interval<T>& w, const _interval<T>& x, const _interval<T>& y)
	// Hull([x], [y])
	{
		T inf, sup;
		inf = ISLIB_MIN(x._inf, y._inf);
		sup = ISLIB_MAX(x._sup, y._sup);
		w = _interval<T>{ inf, sup };
	}

	template <class T>
	void hull(_interval<T>& w, const _interval<T>& x, const T& y)
	// Hull([x], y)
	{
		T inf, sup;
		inf = ISLIB_MIN(x._inf, y);
		sup = ISLIB_MAX(x._sup, y);
		w = _interval<T>{ inf, sup };
	}

	template <class T>
	void hull(_interval<T>& w, const T& x, const T& y)
	// Hull(x, y)
	{
		T inf, sup;
		inf = ISLIB_MIN(x, y);
		sup = ISLIB_MAX(x, y);
		w = _interval<T>{ inf, sup };
	}

	template <class T>
	void hull(_interval<T>& w, const T& y, const _interval<T>& x)
	// Hull(y, [x])
	{
		T inf, sup;
		inf = ISLIB_MIN(x.inf(), y);
		sup = ISLIB_MAX(x.sup(), y);
		w = _interval<T>{ inf, sup };
	}

	template <class T>
	bool _interval<T>::set_contains(const _interval& x) const
	// checks whether an interval [x] is a subinterval of this
	{
		if (x.is_empty()) {
			return true;
		}
		return (x._inf >= _inf) && (x._sup <= _sup);
	}

	template <class T>
	bool _interval<T>::set_contains(const T& x) const
	// checks if a number x is in interval this
	{
		if (_isnan(x)) {
			return false; // true
		}
		return ((x >= _inf) && (x <= _sup)); // zmienione z <= i >=
	}

	template <class T>
	bool _interval<T>::set_strictly_contains(const _interval& x) const
	// checks whether interval [x] a proper subinterval of this
	{
		if (x.is_empty()) {
			return true;
		}
		return ((x._inf > _inf) && (x._sup < _sup));
	}

	template <class T>
	bool _interval<T>::set_strictly_contains(const T& x) const
	// checks if a number x is in interval this
	{
		if (_isnan(x)) {
			return false; //true
		}
		return ((x > _inf) && (x < _sup)); // zmienione z <= i >=
	}

	template <class T>
	// intersection of [x] and [y]
	_interval<T> operator & (const _interval<T>& x, const _interval<T>& y)
	{
		T inf, sup;
		
		if (x.is_empty() || y.is_empty()) {
			return _interval<T>::emptyset();
		}
		inf = ISLIB_MAX(x.inf(), y.inf());
		sup = ISLIB_MIN(x.sup(), y.sup());
		return _interval<T>{ inf, sup };
	}

	template <class T>
	_interval<T> _interval<T>::operator ~ () const
	// dual operator
	{
		if (is_empty()) {
			return (*this);
		}
		return _interval<T>{ _sup, _inf };
	}

	//%/////////////////////////////////////
	//%
	//%  Arithmetic operators
	//%
	//%/////////////////////////////////////

	template <class T>
	// addition [a] + [b]
	_interval<T> operator + (const _interval<T>& x, const _interval<T>& y)
	{
		T inf, sup;
		RoundUp();
		inf = -x.inf() - y.inf();
		sup = x.sup() + y.sup();
		RoundNear();
		return _interval<T>{ -inf, sup };
		//return add(x, y);
	}

	template <class T>
	// [a] + b
	_interval<T> operator + (const _interval<T>& x, const T& y)
	{
		return x + _interval<T>{ y };
		//return add(x, y);
	}

	template <class T>
	// a + [b]
	_interval<T> operator + (const T& y, const _interval<T>& x)
	{
		return _interval<T>{ y } +x;
		//return add(y, x);
	}

	template <class T>
	// [x] += [y]
	void operator += (_interval<T>& x, const _interval<T>& y)
	{
		x = x + y;
		//x = add(x, y);
	}

	template <class T>
	// [x] += y
	void operator += (_interval<T>& x, const T& y)
	{
		x = x + y;
		//x = add(x, y);
	}

	template <class T>
	_interval<T> _interval<T>::operator - () const
	// negation operator
	{
		return _interval<T>{ -_sup, -_inf };
	}

	template <class T>
	// subtraction [x] - [y]
	_interval<T> operator - (const _interval<T>& x, const _interval<T>& y)
	{
		T inf, sup;
		RoundUp();
		inf = -x.inf() + y.sup();
		sup = x.sup() - y.inf();
		RoundNear();
		return _interval<T>{ -inf, sup };
		//return sub(x, y);
	}

	template <class T>
	//% subtraction [x] - y
	_interval<T> operator - (const _interval<T>& x, const T& y)
	{
		return x - _interval<T>{ y };
		//return sub(x, y);
	}

	template <class T>
	//% subtraction x - [y]
	_interval<T> operator - (const T& y, const _interval<T>& x)
	{
		return _interval<T>{ y } -x;
		//return sub(y, x);
	}

	template <class T>
	//% decrementation [x] -= [y]
	void operator -= (_interval<T>& x, const _interval<T>& y)
	{
		x = x - y;
		//x = sub(x, y);
	}

	template <class T>
	//% decrementation [x] -= y
	void operator -= (_interval<T>& x, const T& y)
	{
		x = x - y;
		//x = sub(x, y);
	}

	template <class T>
	//% multiplication of intervals [x] * [y]
	_interval<T> operator * (const _interval<T>& x, const _interval<T>& y)
	{
		T l1, l2, l3, l4, u1, u2, u3, u4;
		T z, z2, z3, z4;
		T xinf = x.inf(), xsup = x.sup(), yinf = y.inf(), ysup = y.sup();

		if (x.is_empty() || y.is_empty()) {
			return _interval<T>::emptyset();
		}
		if (x == 0.0 || y == 0.0) {
			return _interval<T>{ 0.0 };
		}
		if (x.set_contains(0.0) && y.set_contains(0.0)) { // x, y belong to Z = {a\inIR : a.inf * a.sup <= 0}
			l1 = -xinf, l2 = ysup;
			l3 = -xsup, l4 = yinf;
			u1 = xinf, u2 = yinf;
			u3 = xsup, u4 = ysup;
			RoundUp();
			z = l1 * l2;
			z3 = l3 * l4;
			z2 = u1 * u2;
			z4 = u3 * u4;
			RoundNear();
			return _interval<T>{ ISLIB_MIN(-z, -z3), ISLIB_MAX(z2, z4) };
		}
		if (x >= 0.0 && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -yinf, u2 = ysup;
		}
		if (x >= 0.0 && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = ysup;
		}
		if (x <= 0.0 && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -ysup, u2 = yinf;
		}
		if (x <= 0.0 && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -ysup, u2 = yinf;
		}
		if (x >= 0.0 && y.set_contains(0.0)) {
			l1 = xsup, u1 = xsup;
			l2 = -yinf, u2 = ysup;
		}
		if (x.set_contains(0.0) && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -ysup, u2 = ysup;
		}
		if (x <= 0.0 && y.set_contains(0.0)) {
			l1 = xinf, u1 = xinf;
			l2 = -ysup, u2 = yinf;
		}
		if (x.set_contains(0.0) && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = yinf;
		}
		if (x.set_contains(0.0) && y.set_contains(0.0)) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = yinf;
		}
		RoundUp();
		z = l1 * l2;
		z2 = u1 * u2;
		RoundNear();
		return _interval<T>{ -z, z2 }; // required if fstp is used
		//return mult(x, y);
	}

	template <class T>
	//% multiplication [x] * y
	_interval<T> operator * (const _interval<T>& x, const T& y)
	{
		if (x != 0.0) {
			return x * _interval<T>{ y };
		}
		return _interval<T>{ 0.0 };
	}

	template <class T>
	//% multiplication x * [y]
	_interval<T> operator * (const T& y, const _interval<T>& x)
	{
		if (x != 0.0) {
			return _interval<T>{ y } * x;
		}
		return _interval<T>(0.0);
	}

	template <class T>
	_interval<T> operator *= (_interval<T>& x, const _interval<T>& y)
	//% muliplicative incrementation
	{
		if (y == 0.0) {
			x = 0.0;
		}
		else {
			x = x * y;
		}
		return x;
	}

	template <class T>
	_interval<T> operator *= (_interval<T>& x, const T& y)
	// muliplicative incrementation
	{
		if (y == (T) 0.0) {
			x = 0.0;
		}
		else {
			x = x * y;
		}
		return x;
	}

	template <class T>
	_interval<T> _interval<T>::reciprocal() const
	{
		real inf, sup;

		RoundDown();
		inf = 1.0 / _sup;
		RoundUp();
		sup = 1.0 / _inf;
		RoundNear();
		return _interval<T>{ inf, sup };
	}

	template <class T>
	_interval<T> operator / (const _interval<T>& x, const _interval<T>& y)
	{
		if (y.set_contains(0.0)) return _interval<T>::emptyset();
		return x * y.reciprocal();
	}

	template <class T>
	_interval<T> operator / (const _interval<T>& x, const T& y)
	{
		assert(y != 0.0);
		return x * _interval<T>(y).reciprocal();
	}

	template <class T>
	_interval<T> operator / (const T& y, const _interval<T>& x)
	{
		assert(!x.set_contains(0.0));
		return _interval<T>{ y } * x.reciprocal();
	}

	//%//////////////////////////////////////////////////////////////////////
	//%
	//% Comparison operators
	//%
	//%//////////////////////////////////////////////////////////////////////

	template <class T>
	bool operator<(const _interval<T>& x, const _interval<T>& y)
	{
		return (x.sup() < y.inf());
	}

	template <class T>
	bool operator<(const _interval<T>& x, const T& y)
	{
		return (x.sup() < ((_interval<T>) y).inf());
	}

	template <class T>
	bool operator<=(const _interval<T>& x, const _interval<T>& y)
	{
		return (x.inf() <= y.inf() && x.sup() <= y.sup());
	}

	template <class T>
	bool operator<=(const _interval<T>& x, const T& y)
	{
		return (x.sup() <= y);
	}

	template <class T>
	bool operator>(const _interval<T>& x, const _interval<T>& y)
	{
		return (x.inf() > y.sup());
	}

	template <class T>
	bool operator > (const _interval<T>& x, const T& y)
	{
		return (x.inf() > y);
	}

	template <class T>
	bool operator>=(const _interval<T>& x, const _interval<T>& y)
	{
		return (x.inf() >= y.inf() && x.sup() >= y.sup());
	}

	template <class T>
	bool operator>=(const _interval<T>& x, const T& y)
	{
		return (x.inf() >= ((_interval<T>) y).sup());
	}

	template <class T>
	bool operator == (const _interval<T>& x, const _interval<T>& y)
	{
		if (x.is_empty() && y.is_empty()) {
			return true;
		}
		return (x.inf() == y.inf() && x.sup() == y.sup());
	}

	template <class T>
	bool operator == (const _interval<T>& x, const T& y)
	{
		return (x == _interval<T>{ y });
	}

	template <class T>
	bool operator != (const _interval<T>& x, const _interval<T>& y)
	{
		return !(x == y);
	}

	template <class T>
	bool operator != (const _interval<T>& x, const T& y)
	{
		return !(x == y);
	}

	template <class T>
	/* output procedure */
	inline std::ostream& operator << (std::ostream& os, const _interval<T>& x)
	{
		if (x.is_empty()) {
			os << "empty";
			return os;
		}
		os << "[" << x.inf() << ", " << x.sup() << "]";
		return os;
	}

	/////////////////////////////////////////
	//// asm functions
	/////////////////////////////////////////

	/*---------------------------------*
	 * Addition
	 *---------------------------------*/
	template<class T>
	_interval<T> add(const _interval<T>& x,
		const _interval<T>& y)
	{

		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y.inf(), ysup = y.sup();
		short cw;
		short rw;

		__asm {
			fld			xinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		yinf
				fstp		inf
				fld			xsup
				fadd		ysup
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup);
	}

	template<class T>
	_interval<T> add(const _interval<T>& x, const T& y)
	{

		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y, ysup = y;
		short cw;
		short rw;

		__asm {
			fld			xinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		yinf
				fstp		inf
				fld			xsup
				fadd		ysup
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup);
	}

	template<class T>
	_interval<T> add(const T& y, const _interval<T>& x)
	{

		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y, ysup = y;
		short cw;
		short rw;

		__asm {
			fld			yinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		xinf
				fstp		inf
				fld			ysup
				fadd		xsup
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup); // required if fstp is used
	}

	/*---------------------------------*
	 * Subtraction
	 *---------------------------------*/
	template<class T>
	_interval<T> sub(const _interval<T>& x, const _interval<T>& y)
	{
		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y.inf(), ysup = y.sup();
		short cw;
		short rw;

		__asm {
			fld			xinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		ysup
				fstp		inf
				fld			xsup
				fadd		yinf
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup); // required if fstp is used
	}

	template<class T>
	_interval<T> sub(const _interval<T>& x, const T& y)
	{
		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y, ysup = y;
		short cw;
		short rw;

		__asm {
			fld			xinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		ysup
				fstp		inf
				fld			xsup
				fadd		yinf
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup); // required if fstp is used
	}

	template<class T>
	_interval<T> sub(const T& y, const _interval<T>& x)
	{
		T inf, sup;
		T xinf = -x.inf(), xsup = x.sup(), yinf = -y, ysup = y;
		short cw;
		short rw;

		__asm {
			fld			yinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0F3FFh
				or			ax, 00800h
				mov			rw, ax
				fldcw		rw
				fadd		xsup
				fstp		inf
				fld			ysup
				fadd		xinf
				fstp		sup
				fldcw		cw
		}
		return _interval<T>(-inf, sup); // required if fstp is used
	}

	/*---------------------------------*
	 * Multiplication
	 *---------------------------------*/
	template<class T>
	_interval<T> mult(const _interval<T>& x, const _interval<T>& y)
	{
		// case 1: x, y >= 0
		T xinf = x.inf(), yinf = y.inf(), xsup = x.sup(), ysup = y.sup();
		T l1, u1, l2, u2;
		T z, z2;
		short cw;
		short rw;
		//----------------------------
		if (x.set_contains(0.0) && y.set_contains(0.0)) { // x, y belong to Z = {a\inIR : a.inf * a.sup <= 0}
			T l3, l4, u3, u4;
			l1 = -xinf; l2 = ysup;
			l3 = -xsup; l4 = yinf;
			u1 = xinf; u2 = yinf;
			u3 = xsup; u4 = ysup;
			T z3, z4;
			__asm {
				fld			l1 //xinf
					fnstcw		cw
					mov			ax, cw
					and			ax, 0xf3ff //0F3FFh // CRoundNear
					or			ax, 0x0800 //00800h // CRoundUp
					mov			rw, ax
					fldcw		rw
					fmul		l2 //yinf
					fstp		z
					fld			u1 //xsup
					fmul		u2 //ysup
					fstp		z2
					fld			l3
					fmul		l4
					fstp		z3
					fld			u3
					fmul		u4
					fstp		z4
					fldcw		cw
			}
			return _interval<T>(ISLIB_MIN(-z, -z3), ISLIB_MAX(z2, z4));
		}
		if (x >= 0.0 && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -yinf, u2 = ysup;
		}
		if (x >= 0.0 && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = ysup;
		}
		if (x <= 0.0 && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -ysup, u2 = yinf;
		}
		if (x <= 0.0 && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -ysup, u2 = yinf;
		}
		if (x >= 0.0 && y.set_contains(0.0)) {
			l1 = xsup, u1 = xsup;
			l2 = -yinf, u2 = ysup;
		}
		if (x.set_contains(0.0) && y >= 0.0) {
			l1 = xinf, u1 = xsup;
			l2 = -ysup, u2 = ysup;
		}
		if (x <= 0.0 && y.set_contains(0.0)) {
			l1 = xinf, u1 = xinf;
			l2 = -ysup, u2 = yinf;
		}
		if (x.set_contains(0.0) && y <= 0.0) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = yinf;
		}
		if (x.set_contains(0.0) && y.set_contains(0.0)) {
			l1 = xsup, u1 = xinf;
			l2 = -yinf, u2 = yinf;
		}
		//=========================================
		__asm {
			fld			l1 //xinf
				fnstcw		cw
				mov			ax, cw
				and			ax, 0xf3ff //0F3FFh // CRoundNear
				or			ax, 0x0800 //00800h // CRoundUp
				mov			rw, ax
				fldcw		rw
				fmul		l2 //yinf
				fstp		z
				fld			u1 //xsup
				fmul		u2 //ysup
				fstp		z2
				fldcw		cw
		}
		return _interval<T>(-z, z2); // required if fstp is used
	}

	template<class T>
	_interval<T> rcprcl(const _interval<T>& x)
	{
		// reciprocal of an interval
		T xinf = x.inf(), xsup = x.sup(), l = -1.0, l2 = 1.0;
		T z, z2;
		short cw;
		short rw;

		__asm {
			fld			l
				fnstcw		cw
				mov			ax, cw
				and			ax, 0xf3ff //0F3FFh // CRoundNear
				or			ax, 0x0800 //00800h // CRoundUp
				mov			rw, ax
				fldcw		rw
				fdiv		xsup
				fstp		z
				fld			l2 //xsup
				fdiv		xinf //ysup
				fstp		z2
				fldcw		cw
		}
		return _interval<T>(-z, z2); // required if fstp is used
	}

	typedef _interval<real> interval;
}
