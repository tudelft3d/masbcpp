/***********************************************************************
Point - Class for affine points.
Copyright (c) 2001-2010 Oliver Kreylos

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

#ifndef GEOMETRY_POINT_INCLUDED
#define GEOMETRY_POINT_INCLUDED

#include <Math/Math.h>
#include <Geometry/ComponentArray.h>
#include <Geometry/Vector.h>

namespace Geometry {

/* Forward declarations: */
template <class ScalarParam,int dimensionParam>
class AffineCombiner;

template <class ScalarParam,int dimensionParam>
class Point:public Geometry::ComponentArray<ScalarParam,dimensionParam>
	{
	/* Declarations of inherited types/elements: */
	public:
	using ComponentArray<ScalarParam,dimensionParam>::dimension;
	using ComponentArray<ScalarParam,dimensionParam>::components;
	
	/* Embedded classes: */
	public:
	typedef Geometry::Vector<ScalarParam,dimensionParam> Vector; // Compatible vector type
	typedef Geometry::AffineCombiner<ScalarParam,dimensionParam> AffineCombiner; // Compatible affine combiner type
	
	/* Constructors and destructors: */
	static const Point origin; // The zero point (origin of local coordinate system)
	Point(void) // No initialization
		{
		}
	explicit Point(ScalarParam filler) // Fills the component array with a single value
		:ComponentArray<ScalarParam,dimensionParam>(filler)
		{
		}
	Point(ScalarParam c1,ScalarParam c2) // Constructor for 2D point
		:ComponentArray<ScalarParam,dimensionParam>(c1,c2)
		{
		}
	Point(ScalarParam c1,ScalarParam c2,ScalarParam c3) // Constructor for 3D point
		:ComponentArray<ScalarParam,dimensionParam>(c1,c2,c3)
		{
		}
	Point(ScalarParam c1,ScalarParam c2,ScalarParam c3,ScalarParam c4) // Constructor for 4D point
		:ComponentArray<ScalarParam,dimensionParam>(c1,c2,c3,c4)
		{
		}
	template <class SourceScalarParam>
	Point(const SourceScalarParam* array) // Construction from C-style array
		:ComponentArray<ScalarParam,dimensionParam>(array)
		{
		}
	template <class SourceScalarParam,int sourceDimension>
	explicit Point(const ComponentArray<SourceScalarParam,sourceDimension>& source) // Constructs point from component array
		:ComponentArray<ScalarParam,dimensionParam>(source)
		{
		}
	template <class SourceScalarParam,int sourceDimension>
	Point(const Point<SourceScalarParam,sourceDimension>& source) // Copy constructor with type conversion and dimension change
		:ComponentArray<ScalarParam,dimensionParam>(source)
		{
		}
	
	/* Methods: */
	Point operator+(void) const // Unary plus operator; returns copy of point
		{
		return *this;
		}
	Point& operator+=(const Geometry::Vector<ScalarParam,dimensionParam>& v) // Addition assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]+=v[i];
		return *this;
		}
	Point& operator-=(const Geometry::Vector<ScalarParam,dimensionParam>& v) // Subtraction assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]-=v[i];
		return *this;
		}
	};

/************************************
Operations on objects of class Point:
************************************/

template <class ScalarParam,int dimensionParam>
inline Point<ScalarParam,dimensionParam> operator+(const Point<ScalarParam,dimensionParam>& p,const Vector<ScalarParam,dimensionParam>& v) // Addition of point and vector
	{
	Point<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=p[i]+v[i];
	return result;
	}

template <class ScalarParam>
inline Point<ScalarParam,2> operator+(const Point<ScalarParam,2>& p,const Vector<ScalarParam,2>& v)
	{
	return Point<ScalarParam,2>(p[0]+v[0],p[1]+v[1]);
	}

template <class ScalarParam>
inline Point<ScalarParam,3> operator+(const Point<ScalarParam,3>& p,const Vector<ScalarParam,3>& v)
	{
	return Point<ScalarParam,3>(p[0]+v[0],p[1]+v[1],p[2]+v[2]);
	}

template <class ScalarParam,int dimensionParam>
inline Point<ScalarParam,dimensionParam> operator+(const Vector<ScalarParam,dimensionParam>& v,const Point<ScalarParam,dimensionParam>& p) // Ditto
	{
	Point<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=v[i]+p[i];
	return result;
	}

template <class ScalarParam>
inline Point<ScalarParam,2> operator+(const Vector<ScalarParam,2>& v,const Point<ScalarParam,2>& p)
	{
	return Point<ScalarParam,2>(v[0]+p[0],v[1]+p[1]);
	}

template <class ScalarParam>
inline Point<ScalarParam,3> operator+(const Vector<ScalarParam,3>& v,const Point<ScalarParam,3>& p)
	{
	return Point<ScalarParam,3>(v[0]+p[0],v[1]+p[1],v[2]+p[2]);
	}

template <class ScalarParam,int dimensionParam>
inline Point<ScalarParam,dimensionParam> operator-(const Point<ScalarParam,dimensionParam>& p,const Vector<ScalarParam,dimensionParam>& v) // Subtraction of point and vector
	{
	Point<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=p[i]-v[i];
	return result;
	}

template <class ScalarParam>
inline Point<ScalarParam,2> operator-(const Point<ScalarParam,2>& p,const Vector<ScalarParam,2>& v)
	{
	return Point<ScalarParam,2>(p[0]-v[0],p[1]-v[1]);
	}

template <class ScalarParam>
inline Point<ScalarParam,3> operator-(const Point<ScalarParam,3>& p,const Vector<ScalarParam,3>& v)
	{
	return Point<ScalarParam,3>(p[0]-v[0],p[1]-v[1],p[2]-v[2]);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator-(const Point<ScalarParam,dimensionParam>& p1,const Point<ScalarParam,dimensionParam>& p2) // Distance vector between two points
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=p1[i]-p2[i];
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator-(const Point<ScalarParam,2>& p1,const Point<ScalarParam,2>& p2)
	{
	return Vector<ScalarParam,2>(p1[0]-p2[0],p1[1]-p2[1]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator-(const Point<ScalarParam,3>& p1,const Point<ScalarParam,3>& p2)
	{
	return Vector<ScalarParam,3>(p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]);
	}

template <class ScalarParam,int dimensionParam>
inline ScalarParam sqrDist(const Point<ScalarParam,dimensionParam>& p1,const Point<ScalarParam,dimensionParam>& p2) // Squared Euklidean distance between two points
	{
	ScalarParam result(0);
	for(int i=0;i<dimensionParam;++i)
		result+=Math::sqr(p1[i]-p2[i]);
	return result;
	}

template <class ScalarParam>
inline ScalarParam sqrDist(const Point<ScalarParam,2>& p1,const Point<ScalarParam,2>& p2)
	{
	return Math::sqr(p1[0]-p2[0])+Math::sqr(p1[1]-p2[1]);
	}

template <class ScalarParam>
inline ScalarParam sqrDist(const Point<ScalarParam,3>& p1,const Point<ScalarParam,3>& p2)
	{
	return Math::sqr(p1[0]-p2[0])+Math::sqr(p1[1]-p2[1])+Math::sqr(p1[2]-p2[2]);
	}

template <class ScalarParam,int dimensionParam>
inline double dist(const Point<ScalarParam,dimensionParam>& p1,const Point<ScalarParam,dimensionParam>& p2) // Euklidean distance between two points
	{
	double result=0.0;
	for(int i=0;i<dimensionParam;++i)
		result+=Math::sqr(double(p1[i])-double(p2[i]));
	return Math::sqrt(result);
	}

template <class ScalarParam>
inline double dist(const Point<ScalarParam,2>& p1,const Point<ScalarParam,2>& p2)
	{
	return Math::sqrt(Math::sqr(double(p1[0])-double(p2[0]))+Math::sqr(double(p1[1])-double(p2[1])));
	}

template <class ScalarParam>
inline double dist(const Point<ScalarParam,3>& p1,const Point<ScalarParam,3>& p2)
	{
	return Math::sqrt(Math::sqr(double(p1[0])-double(p2[0]))+Math::sqr(double(p1[1])-double(p2[1]))+Math::sqr(double(p1[2])-double(p2[2])));
	}

template <class ScalarParam,int dimensionParam>
inline Point<ScalarParam,dimensionParam> mid(const Point<ScalarParam,dimensionParam>& p1,const Point<ScalarParam,dimensionParam>& p2) // Returns midpoint of two points
	{
	Point<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=Math::mid(p1[i],p2[i]);
	return result;
	}

template <class ScalarParam>
inline Point<ScalarParam,2> mid(const Point<ScalarParam,2>& p1,const Point<ScalarParam,2>& p2)
	{
	return Point<ScalarParam,2>(Math::mid(p1[0],p2[0]),Math::mid(p1[1],p2[1]));
	}

template <class ScalarParam>
inline Point<ScalarParam,3> mid(const Point<ScalarParam,3>& p1,const Point<ScalarParam,3>& p2)
	{
	return Point<ScalarParam,3>(Math::mid(p1[0],p2[0]),Math::mid(p1[1],p2[1]),Math::mid(p1[2],p2[2]));
	}

template <class ScalarParam,int dimensionParam>
inline Point<ScalarParam,dimensionParam> affineCombination(const Point<ScalarParam,dimensionParam>& p1,const Point<ScalarParam,dimensionParam>& p2,ScalarParam w2) // Affine combination of two points
	{
	ScalarParam w1=ScalarParam(1)-w2;
	Point<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=p1[i]*w1+p2[i]*w2;
	return result;
	}

template <class ScalarParam>
inline Point<ScalarParam,2> affineCombination(const Point<ScalarParam,2>& p1,const Point<ScalarParam,2>& p2,ScalarParam w2)
	{
	ScalarParam w1=ScalarParam(1)-w2;
	return Point<ScalarParam,2>(p1[0]*w1+p2[0]*w2,p1[1]*w1+p2[1]*w2);
	}

template <class ScalarParam>
inline Point<ScalarParam,3> affineCombination(const Point<ScalarParam,3>& p1,const Point<ScalarParam,3>& p2,ScalarParam w2)
	{
	ScalarParam w1=ScalarParam(1)-w2;
	return Point<ScalarParam,3>(p1[0]*w1+p2[0]*w2,p1[1]*w1+p2[1]*w2,p1[2]*w1+p2[2]*w2);
	}

}

#if defined(GEOMETRY_NONSTANDARD_TEMPLATES) && !defined(GEOMETRY_POINT_IMPLEMENTATION)
#include <Geometry/Point.icpp>
#endif

#endif
