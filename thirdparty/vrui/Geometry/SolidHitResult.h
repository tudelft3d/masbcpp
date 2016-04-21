/***********************************************************************
SolidHitResult - Class to report results of intersection tests between
rays and solid geometric objects.
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

#ifndef GEOMETRY_SOLIDHITRESULT_INCLUDED
#define GEOMETRY_SOLIDHITRESULT_INCLUDED

#include <Geometry/HitResult.h>

namespace Geometry {

template <class ScalarParam>
class SolidHitResult:public HitResult<ScalarParam> // Class to report intersections with closed surfaces of solids
	{
	/* Embedded classes: */
	public:
	enum Direction // Intersection "direction,", i.e., entering or leaving the solid
		{
		INVALID_DIRECTION,ENTRY,EXIT
		};
	
	/* Elements: */
	private:
	Direction direction;
	
	/* Constructors and destructors: */
	public:
	SolidHitResult(void) // Constructs invalid hit result
		:direction(INVALID_DIRECTION)
		{
		}
	SolidHitResult(ScalarParam sLambda,Direction sDirection) // Constructs valid hit result
		:HitResult<ScalarParam>(sLambda),direction(sDirection)
		{
		}
	
	/* Methods: */
	Direction getDirection(void) const // Returns direction of reported intersection
		{
		return direction;
		}
	};

}

#endif
