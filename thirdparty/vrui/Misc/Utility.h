/***********************************************************************
Utility - Helper functions and classes for a variety of typical tasks.
Copyright (c) 2007 Oliver Kreylos

This file is part of the Miscellaneous Support Library (Misc).

The Miscellaneous Support Library is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

The Miscellaneous Support Library is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Miscellaneous Support Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#ifndef MISC_UTILITY_INCLUDED
#define MISC_UTILITY_INCLUDED

namespace Misc {

/**********************************
Swap two values of arbitrary types:
**********************************/

template <class ValueParam>
inline
void
swap(
	ValueParam& v1,
	ValueParam& v2)
	{
	ValueParam temp=v1;
	v1=v2;
	v2=temp;
	}

/**************************************************************
Calculate the minimum or maximum of two values with operator<=:
**************************************************************/

template <class ValueParam>
inline
const ValueParam&
min(
	const ValueParam& v1,
	const ValueParam& v2)
	{
	if(v1<=v2)
		return v1;
	else
		return v2;
	}

template <class ValueParam>
inline
const ValueParam&
max(
	const ValueParam& v1,
	const ValueParam& v2)
	{
	if(v1<=v2)
		return v2;
	else
		return v1;
	}

}

#endif
