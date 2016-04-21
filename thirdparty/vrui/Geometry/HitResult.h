/***********************************************************************
HitResult - Base classes to report results of intersection tests between
rays and geometric objects. Geometry object classes are expected to
overload this class to return additional information specific to them.
Copyright (c) 2002-2005 Oliver Kreylos

This file is part of the Templatized Geometry Library (TGL).

The Templatized Geometry Library is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The Templatized Geometry Library is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Templatized Geometry Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#ifndef GEOMETRY_HITRESULT_INCLUDED
#define GEOMETRY_HITRESULT_INCLUDED

#include <Math/Constants.h>

namespace Geometry {

template <class ScalarParam>
class HitResult // Base class for intersections between rays and general surfaces
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // Type for scalar values (intersection parameters)
	
	/* Elements: */
	private:
	bool hit; // Flag if intersection is valid
	Scalar lambda; // Ray parameter of reported intersection point
	
	/* Constructors and destructors: */
	public:
	HitResult(void) // Constructs invalid hit result
		:hit(false),lambda(Math::Constants<Scalar>::max)
		{
		}
	HitResult(Scalar sLambda) // Constructs valid hit result
		:hit(true),lambda(sLambda)
		{
		}
	
	/* Methods: */
	bool isValid(void) const // Checks if hit result is valid
		{
		return hit;
		}
	Scalar getParameter(void) const // Returns ray parameter of intersection point
		{
		return lambda;
		}
	};

}

#endif
