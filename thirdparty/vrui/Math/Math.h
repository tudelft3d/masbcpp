/***********************************************************************
Math - Genericized versions of standard C math functions.
Copyright (c) 2001-2012 Oliver Kreylos

This file is part of the Templatized Math Library (Math).

The Templatized Math Library is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The Templatized Math Library is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Templatized Math Library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
***********************************************************************/

#ifndef MATH_MATH_INCLUDED
#define MATH_MATH_INCLUDED

#include <math.h>
#include <stdlib.h>

/* Check if the implementation provides float versions of math calls: */
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
#define MATH_CONFIG_HAVE_FLOAT_CALLS
#endif

/* Check if the implementation provides float classification functions: */
#if defined(__GNUC__) && !defined(__APPLE__)
#define MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
#endif

namespace Math {

/**********************************************
Floating-point number classification functions:
**********************************************/

template <class ScalarParam>
inline bool isNan(ScalarParam value)
	{
	/* General types don't have NAN: */
	return false;
	}

template <>
inline bool isNan(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __isnanf(value);
	#else
	return isnan(value);
	#endif
	}

template <>
inline bool isNan(double value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __isnan(value);
	#else
	return isnan(value);
	#endif
	}

template <class ScalarParam>
inline bool isInf(ScalarParam value)
	{
	/* General types don't have infinity: */
	return false;
	}

template <>
inline bool isInf(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __isinff(value);
	#else
	return isinf(value);
	#endif
	}

template <>
inline bool isInf(double value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __isinf(value);
	#else
	return isinf(value);
	#endif
	}

template <class ScalarParam>
inline bool isFinite(ScalarParam value)
	{
	/* General types are always finite: */
	return true;
	}

template <>
inline bool isFinite(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __finitef(value);
	#else
	return isfinite(value);
	#endif
	}

template <>
inline bool isFinite(double value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CLASSIFICATIONS
	return __finite(value);
	#else
	return isfinite(value);
	#endif
	}

/******************************************
Optimized arithmetic convenience functions:
******************************************/

/* The copysign function returns a value that has the absolute value of abs, but the sign of sign: */

inline signed char copysign(signed char abs,signed char sign)
	{
	return ((abs^sign)&0x80)?-abs:abs;
	}

inline short copysign(short abs,short sign)
	{
	return ((abs^sign)&0x8000)?-abs:abs;
	}

inline int copysign(int abs,int sign)
	{
	return ((abs^sign)&0x80000000)?-abs:abs;
	}

inline long copysign(long abs,long sign)
	{
	return ((abs^sign)&0x8000000000000000L)?-abs:abs;
	}

inline float copysign(float abs,float sign)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return copysignf(abs,sign);
	#else
	return float(::copysign(double(abs),double(sign)));
	#endif
	}

inline double copysign(double abs,double sign)
	{
	return ::copysign(abs,sign);
	}

template <class ScalarParam>
inline ScalarParam mul2(ScalarParam value)
	{
	return value+value;
	}

template <class ScalarParam>
inline ScalarParam div2(ScalarParam value)
	{
	return value/ScalarParam(2);
	}

template <>
inline float div2(float value)
	{
	return value*0.5f;
	}

template <>
inline double div2(double value)
	{
	return value*0.5;
	}

template <class ScalarParam>
inline ScalarParam mid(ScalarParam value1,ScalarParam value2)
	{
	return (value1+value2)/ScalarParam(2);
	}

inline float mid(float value1,float value2)
	{
	return (value1+value2)*0.5f;
	}

inline double mid(double value1,double value2)
	{
	return (value1+value2)*0.5;
	}

template <class ScalarParam>
inline ScalarParam sqr(ScalarParam value)
	{
	return value*value;
	}

template <class ScalarParam>
inline ScalarParam min(ScalarParam v1,ScalarParam v2)
	{
	return v1<=v2?v1:v2;
	}

template <class ScalarParam>
inline ScalarParam max(ScalarParam v1,ScalarParam v2)
	{
	return v1>=v2?v1:v2;
	}

template <class ScalarParam>
inline ScalarParam clamp(ScalarParam value,ScalarParam min,ScalarParam max)
	{
	/* Limit the value to the valid range: */
	if(value<min)
		value=min;
	if(value>max)
		value=max;
	
	/* Return the potentially modified value: */
	return value;
	}

/*************************************************
Type-safe wrappers around standard math functions:
*************************************************/

inline int abs(int value)
	{
	return ::abs(value);
	}

inline float abs(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return fabsf(value);
	#else
	return float(fabs(double(value)));
	#endif
	}

inline double abs(double value)
	{
	return fabs(value);
	}

inline int floor(int value)
	{
	return value;
	}

inline float floor(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return floorf(value);
	#else
	return float(::floor(double(value)));
	#endif
	}

inline double floor(double value)
	{
	return ::floor(value);
	}

inline int ceil(int value)
	{
	return value;
	}

inline float ceil(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return ceilf(value);
	#else
	return float(::ceil(double(value)));
	#endif
	}

inline double ceil(double value)
	{
	return ::ceil(value);
	}

inline int mod(int counter,int denominator)
	{
	return counter%denominator;
	}

inline float mod(float counter,float denominator)
	{
	return float(fmod(counter,denominator));
	}

inline double mod(double counter,double denominator)
	{
	return fmod(counter,denominator);
	}

template <class ScalarParam>
inline ScalarParam rem(ScalarParam counter,ScalarParam denominator)
	{
	ScalarParam result=mod(counter,denominator);
	if(result<ScalarParam(0))
		result+=denominator;
	return result;
	}

inline float sqrt(float value)
	{
	return float(::sqrt(double(value)));
	}

inline double sqrt(double value)
	{
	return ::sqrt(value);
	}

/*********************************
Helper functions for trigonometry:
*********************************/

inline float deg(float radians)
	{
	return radians*(180.0f/3.14159265358979323846f);
	}

inline double deg(double radians)
	{
	return radians*(180.0/3.14159265358979323846);
	}

inline float rad(float degrees)
	{
	return degrees*(3.14159265358979323846f/180.0f);
	}

inline double rad(double degrees)
	{
	return degrees*(3.14159265358979323846/180.0);
	}

inline float wrapRad(float radians)
	{
	return radians-float(floor(double(radians)/(2.0*3.14159265358979323846))*(2.0*3.14159265358979323846));
	}

inline double wrapRad(double radians)
	{
	return radians-floor(radians/(2.0*3.14159265358979323846))*(2.0*3.14159265358979323846);
	}

/*************************************************
Type-safe wrappers around trigonometric functions:
*************************************************/

inline float sin(float radians)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return sinf(radians);
	#else
	return float(::sin(double(radians)));
	#endif
	}

inline double sin(double radians)
	{
	return ::sin(radians);
	}

inline float cos(float radians)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return cosf(radians);
	#else
	return float(::cos(double(radians)));
	#endif
	}

inline double cos(double radians)
	{
	return ::cos(radians);
	}

inline float tan(float radians)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return tanf(radians);
	#else
	return float(::tan(double(radians)));
	#endif
	}

inline double tan(double radians)
	{
	return ::tan(radians);
	}

inline float asin(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return asinf(value);
	#else
	return float(::asin(double(value)));
	#endif
	}

inline double asin(double value)
	{
	return ::asin(value);
	}

inline float acos(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return acosf(value);
	#else
	return float(::acos(double(value)));
	#endif
	}

inline double acos(double value)
	{
	return ::acos(value);
	}

inline float atan(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return atanf(value);
	#else
	return float(::atan(double(value)));
	#endif
	}

inline double atan(double value)
	{
	return ::atan(value);
	}

inline float atan2(float counter,float denominator)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return atan2f(counter,denominator);
	#else
	return float(::atan2(double(counter),double(denominator)));
	#endif
	}

inline double atan2(double counter,double denominator)
	{
	return ::atan2(counter,denominator);
	}

/************************************************************
Type-safe wrappers around hyperbolic trigonometric functions:
************************************************************/

inline float sinh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return sinhf(value);
	#else
	return float(::sinh(double(value)));
	#endif
	}

inline double sinh(double value)
	{
	return ::sinh(value);
	}

inline float cosh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return coshf(value);
	#else
	return float(::cosh(double(value)));
	#endif
	}

inline double cosh(double value)
	{
	return ::cosh(value);
	}

inline float tanh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return tanhf(value);
	#else
	return float(::tanh(double(value)));
	#endif
	}

inline double tanh(double value)
	{
	return ::tanh(value);
	}

inline float asinh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return asinhf(value);
	#else
	return float(::asinh(double(value)));
	#endif
	}

inline double asinh(double value)
	{
	return ::asinh(value);
	}

inline float acosh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return acoshf(value);
	#else
	return float(::acosh(double(value)));
	#endif
	}

inline double acosh(double value)
	{
	return ::acosh(value);
	}

inline float atanh(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return atanhf(value);
	#else
	return float(::atanh(double(value)));
	#endif
	}

inline double atanh(double value)
	{
	return ::atanh(value);
	}

/**************************************************
Type-safe wrappers around transcendental functions:
**************************************************/

inline float log(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return logf(value);
	#else
	return float(::log(double(value)));
	#endif
	}

inline double log(double value)
	{
	return ::log(value);
	}

inline float log10(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return log10f(value);
	#else
	return float(::log10(double(value)));
	#endif
	}

inline double log10(double value)
	{
	return ::log10(value);
	}

inline float exp(float value)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return expf(value);
	#else
	return float(::exp(double(value)));
	#endif
	}

inline double exp(double value)
	{
	return ::exp(value);
	}

inline float pow(float base,float exponent)
	{
	#ifdef MATH_CONFIG_HAVE_FLOAT_CALLS
	return powf(base,exponent);
	#else
	return float(::pow(double(base),double(exponent)));
	#endif
	}

inline double pow(double base,double exponent)
	{
	return ::pow(base,exponent);
	}

}

#endif
