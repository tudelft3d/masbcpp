/***********************************************************************
Ray - Class for affine rays.
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

#ifndef GEOMETRY_RAY_INCLUDED
#define GEOMETRY_RAY_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>

namespace Geometry {

template <class ScalarParam,int dimensionParam>
class Ray
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=dimensionParam; // The ray's dimension
	typedef Geometry::Vector<ScalarParam,dimensionParam> Vector; // The type for vectors
	typedef Geometry::Point<ScalarParam,dimensionParam> Point; // The type for points
	
	/* Elements: */
	private:
	Point origin; // Origin point of ray
	Vector direction; // Direction vector of ray (not normalized)
	
	/* Constructors and destructors: */
	public:
	Ray(void) // Dummy constructor
		{
		}
	Ray(const Point& sOrigin,const Vector& sDirection) // Elementwise
		:origin(sOrigin),direction(sDirection)
		{
		}
	Ray(const Point& p1,const Point& p2) // Ray from point p1 to point p2
		:origin(p1),direction(p2-p1)
		{
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	Ray(const Ray<SourceScalarParam,sourceDimensionParam>& source) // Copy constructor with type conversion
		:origin(source.getOrigin()),direction(source.getDirection())
		{
		}
	
	/* Methods: */
	const Point& getOrigin(void) const // Returns ray's origin
		{
		return origin;
		}
	Ray& setOrigin(const Point& newOrigin) // Sets a new origin
		{
		origin=newOrigin;
		return *this;
		}
	const Vector& getDirection(void) const // Returns ray's direction
		{
		return direction;
		}
	Ray& setDirection(const Vector& newDirection) // Sets a new direction
		{
		direction=newDirection;
		return *this;
		}
	Point operator()(Scalar lambda) const // Evaluates ray for parameter lambda
		{
		return origin+direction*lambda;
		}
	double getDirectionMag(void) const // Returns length of direction vector
		{
		return direction.mag();
		}
	Ray& normalizeDirection(void) // Normalizes direction vector
		{
		direction.normalize();
		return *this;
		}
	template <class TransformationParam>
	Ray& transform(const TransformationParam& t) // Transforms ray
		{
		origin=t.transform(origin);
		direction=t.transform(direction);
		return *this;
		}
	template <class TransformationParam>
	Ray& inverseTransform(const TransformationParam& t) // Transforms ray by the inverse of the given transformation
		{
		origin=t.inverseTransform(origin);
		direction=t.inverseTransform(direction);
		return *this;
		}
	};

}

#endif
