/***********************************************************************
Constants - Classes providing generic access to type-specific math-
relevant information.
Copyright (c) 2003-2005 Oliver Kreylos

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

#ifndef MATH_CONSTANTS_INCLUDED
#define MATH_CONSTANTS_INCLUDED

#ifdef __APPLE__
#include <machine/limits.h>
#endif
#ifdef __linux__
#include <limits.h>
#endif

namespace Math {

/************************************************************
"Dummy" generic class - only specialized versions make sense:
************************************************************/

template <class ScalarParam>
class Constants
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar;
	
	/* Elements: */
	static const bool isIntegral=false; // Not an integer type
	static const bool isRing=false; // Not even an approximation of a ring
	static const bool isField=false; // Not an approximation of a field either
	static const bool isReal=false; // Not a model of real numbers
	};

/***********************************************
Specialized constants for all atomic data types:
***********************************************/

template <>
class Constants<bool>
	{
	/* Embedded classes: */
	public:
	typedef bool Scalar;
	
	/* Elements: */
	static const bool isIntegral=false;
	static const bool isRing=false;
	static const bool isField=false;
	static const bool isReal=false;
	};

template <>
class Constants<signed char>
	{
	/* Embedded classes: */
	public:
	typedef signed char Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=true;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=(signed char)0;
	static const Scalar one=(signed char)1;
	static const Scalar min=SCHAR_MIN;
	static const Scalar max=SCHAR_MAX;
	};

template <>
class Constants<unsigned char>
	{
	/* Embedded classes: */
	public:
	typedef unsigned char Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=false;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=(unsigned char)0;
	static const Scalar one=(unsigned char)1;
	static const Scalar min=(unsigned char)0;
	static const Scalar max=UCHAR_MAX;
	};

template <>
class Constants<char>
	{
	/* Embedded classes: */
	public:
	typedef char Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=false; // This would be true if char==signed char. But who knows?
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=(char)0;
	static const Scalar one=(char)1;
	static const Scalar min=CHAR_MIN;
	static const Scalar max=CHAR_MAX;
	};

template <>
class Constants<short>
	{
	/* Embedded classes: */
	public:
	typedef short Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=true;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=(short)0;
	static const Scalar one=(short)1;
	static const Scalar min=SHRT_MIN;
	static const Scalar max=SHRT_MAX;
	};

template <>
class Constants<unsigned short>
	{
	/* Embedded classes: */
	public:
	typedef unsigned short Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=false;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=(unsigned short)0;
	static const Scalar one=(unsigned short)1;
	static const Scalar min=(unsigned short)0;
	static const Scalar max=USHRT_MAX;
	};

template <>
class Constants<int>
	{
	/* Embedded classes: */
	public:
	typedef int Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=true;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=0;
	static const Scalar one=1;
	static const Scalar min=INT_MIN;
	static const Scalar max=INT_MAX;
	};

template <>
class Constants<unsigned int>
	{
	/* Embedded classes: */
	public:
	typedef unsigned int Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=false;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=0U;
	static const Scalar one=1U;
	static const Scalar min=0U;
	static const Scalar max=UINT_MAX;
	};

template <>
class Constants<long>
	{
	/* Embedded classes: */
	public:
	typedef long Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=true;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=0L;
	static const Scalar one=1L;
	static const Scalar min=LONG_MIN;
	static const Scalar max=LONG_MAX;
	};

template <>
class Constants<unsigned long>
	{
	/* Embedded classes: */
	public:
	typedef unsigned long Scalar;
	
	/* Elements: */
	static const bool isIntegral=true;
	static const bool isRing=false;
	static const bool isField=false;
	static const bool isReal=false;
	static const Scalar zero=0UL;
	static const Scalar one=1UL;
	static const Scalar min=0UL;
	static const Scalar max=ULONG_MAX;
	};

template <>
class Constants<float>
	{
	/* Embedded classes: */
	public:
	typedef float Scalar;
	
	/* Elements: */
	static const bool isIntegral=false;
	static const bool isRing=true;
	static const bool isField=true;
	static const bool isReal=true;
	static const Scalar zero;
	static const Scalar one;
	static const Scalar min;
	static const Scalar max;
	static const Scalar smallest;
	static const Scalar epsilon;
	static const Scalar e;
	static const Scalar pi;
	};

template <>
class Constants<double>
	{
	/* Embedded classes: */
	public:
	typedef double Scalar;
	
	/* Elements: */
	static const bool isIntegral=false;
	static const bool isRing=true;
	static const bool isField=true;
	static const bool isReal=true;
	static const Scalar zero;
	static const Scalar one;
	static const Scalar min;
	static const Scalar max;
	static const Scalar smallest;
	static const Scalar epsilon;
	static const Scalar e;
	static const Scalar pi;
	};

}

#endif
