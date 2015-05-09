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

#include <float.h>

#include <Math/Constants.h>

namespace Math {

/***************************************
Static members of class Constants<bool>:
***************************************/

const bool Constants<bool>::isIntegral;
const bool Constants<bool>::isRing;
const bool Constants<bool>::isField;
const bool Constants<bool>::isReal;

/**********************************************
Static members of class Constants<signed char>:
**********************************************/

const bool Constants<signed char>::isIntegral;
const bool Constants<signed char>::isRing;
const bool Constants<signed char>::isField;
const bool Constants<signed char>::isReal;
const signed char Constants<signed char>::zero;
const signed char Constants<signed char>::one;
const signed char Constants<signed char>::min;
const signed char Constants<signed char>::max;

/************************************************
Static members of class Constants<unsigned char>:
************************************************/

const bool Constants<unsigned char>::isIntegral;
const bool Constants<unsigned char>::isRing;
const bool Constants<unsigned char>::isField;
const bool Constants<unsigned char>::isReal;
const unsigned char Constants<unsigned char>::zero;
const unsigned char Constants<unsigned char>::one;
const unsigned char Constants<unsigned char>::min;
const unsigned char Constants<unsigned char>::max;

/***************************************
Static members of class Constants<char>:
***************************************/

const bool Constants<char>::isIntegral;
const bool Constants<char>::isRing;
const bool Constants<char>::isField;
const bool Constants<char>::isReal;
const char Constants<char>::zero;
const char Constants<char>::one;
const char Constants<char>::min;
const char Constants<char>::max;

/****************************************
Static members of class Constants<short>:
****************************************/

const bool Constants<short>::isIntegral;
const bool Constants<short>::isRing;
const bool Constants<short>::isField;
const bool Constants<short>::isReal;
const short Constants<short>::zero;
const short Constants<short>::one;
const short Constants<short>::min;
const short Constants<short>::max;

/*************************************************
Static members of class Constants<unsigned short>:
*************************************************/

const bool Constants<unsigned short>::isIntegral;
const bool Constants<unsigned short>::isRing;
const bool Constants<unsigned short>::isField;
const bool Constants<unsigned short>::isReal;
const unsigned short Constants<unsigned short>::zero;
const unsigned short Constants<unsigned short>::one;
const unsigned short Constants<unsigned short>::min;
const unsigned short Constants<unsigned short>::max;

/**************************************
Static members of class Constants<int>:
**************************************/

const bool Constants<int>::isIntegral;
const bool Constants<int>::isRing;
const bool Constants<int>::isField;
const bool Constants<int>::isReal;
const int Constants<int>::zero;
const int Constants<int>::one;
const int Constants<int>::min;
const int Constants<int>::max;

/***********************************************
Static members of class Constants<unsigned int>:
***********************************************/

const bool Constants<unsigned int>::isIntegral;
const bool Constants<unsigned int>::isRing;
const bool Constants<unsigned int>::isField;
const bool Constants<unsigned int>::isReal;
const unsigned int Constants<unsigned int>::zero;
const unsigned int Constants<unsigned int>::one;
const unsigned int Constants<unsigned int>::min;
const unsigned int Constants<unsigned int>::max;

/***************************************
Static members of class Constants<long>:
***************************************/

const bool Constants<long>::isIntegral;
const bool Constants<long>::isRing;
const bool Constants<long>::isField;
const bool Constants<long>::isReal;
const long Constants<long>::zero;
const long Constants<long>::one;
const long Constants<long>::min;
const long Constants<long>::max;

/************************************************
Static members of class Constants<unsigned long>:
************************************************/

const bool Constants<unsigned long>::isIntegral;
const bool Constants<unsigned long>::isRing;
const bool Constants<unsigned long>::isField;
const bool Constants<unsigned long>::isReal;
const unsigned long Constants<unsigned long>::zero;
const unsigned long Constants<unsigned long>::one;
const unsigned long Constants<unsigned long>::min;
const unsigned long Constants<unsigned long>::max;

/****************************************
Static members of class Constants<float>:
****************************************/

const bool Constants<float>::isIntegral;
const bool Constants<float>::isRing;
const bool Constants<float>::isField;
const bool Constants<float>::isReal;
const Constants<float>::Scalar Constants<float>::zero=0.0f;
const Constants<float>::Scalar Constants<float>::one=1.0f;
const Constants<float>::Scalar Constants<float>::min=-FLT_MAX;
const Constants<float>::Scalar Constants<float>::max=FLT_MAX;
const Constants<float>::Scalar Constants<float>::smallest=FLT_MIN;
const Constants<float>::Scalar Constants<float>::epsilon=FLT_EPSILON;
const Constants<float>::Scalar Constants<float>::e=2.7182818284590452354f;
const Constants<float>::Scalar Constants<float>::pi=3.14159265358979323846f;

/*****************************************
Static members of class Constants<double>:
*****************************************/

const bool Constants<double>::isIntegral;
const bool Constants<double>::isRing;
const bool Constants<double>::isField;
const bool Constants<double>::isReal;
const Constants<double>::Scalar Constants<double>::zero=0.0;
const Constants<double>::Scalar Constants<double>::one=1.0;
const Constants<double>::Scalar Constants<double>::min=-DBL_MAX;
const Constants<double>::Scalar Constants<double>::max=DBL_MAX;
const Constants<double>::Scalar Constants<double>::smallest=DBL_MIN;
const Constants<double>::Scalar Constants<double>::epsilon=DBL_EPSILON;
const Constants<double>::Scalar Constants<double>::e=2.7182818284590452354;
const Constants<double>::Scalar Constants<double>::pi=3.14159265358979323846;

}
