/***********************************************************************
Box - Class for n-dimensional axis-aligned boxes.
Copyright (c) 2001-2013 Oliver Kreylos

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

#ifndef GEOMETRY_BOX_INCLUDED
#define GEOMETRY_BOX_INCLUDED

#include <utility>
#include <Geometry/ComponentArray.h>
#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Ray.h>
#include <Geometry/SolidHitResult.h>

namespace Geometry {

/* Forward declarations for friend functions: */
template <class ScalarParam,int dimensionParam>
class Box;
template <class ScalarParam,int dimensionParam>
bool operator==(const Box<ScalarParam,dimensionParam>&,const Box<ScalarParam,dimensionParam>&);
template <class ScalarParam,int dimensionParam>
bool operator!=(const Box<ScalarParam,dimensionParam>&,const Box<ScalarParam,dimensionParam>&);
template <class ScalarParam,int dimensionParam>
Box<ScalarParam,dimensionParam> add(const Box<ScalarParam,dimensionParam>&,const Box<ScalarParam,dimensionParam>&);
template <class ScalarParam,int dimensionParam>
Box<ScalarParam,dimensionParam> intersect(const Box<ScalarParam,dimensionParam>&,const Box<ScalarParam,dimensionParam>&);

template <class ScalarParam,int dimensionParam>
class Box
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=dimensionParam; // The box's dimension
	typedef Geometry::Point<ScalarParam,dimensionParam> Point; // The type for points
	typedef Geometry::Vector<ScalarParam,dimensionParam> Vector; // The type for vectors
	typedef Geometry::ComponentArray<ScalarParam,dimensionParam> Size; // The type for size values
	typedef Geometry::Ray<ScalarParam,dimensionParam> Ray; // Compatible ray type
	typedef SolidHitResult<ScalarParam> HitResult; // Hit result type
	
	/* Elements: */
	public:
	Point min,max; // Minimum and maximum point along all primary axes
	
	/* Constructors and destructors: */
	public:
	Box(ScalarParam sMin,ScalarParam sMax) // Used to efficiently construct "special" boxes
		:min(sMin),max(sMax)
		{
		}
	static const Box empty; // The box containing no points at all
	static const Box full; // The box containing every point
	Box(void) // Creates uninitialized box
		{
		}
	Box(const Point& sMin,const Point& sMax) // Constructs a box from a pair of minimal and maximal points
		:min(sMin),max(sMax)
		{
		}
	Box(const Point& sOrigin,const Size& sSize); // Constructs a box from origin and size
	template <class SourceScalarParam,int sourceDimensionParam>
	Box(const Box<SourceScalarParam,sourceDimensionParam>& source) // Copy constructor with type conversion and dimension change
		:min(source.min),max(source.max)
		{
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	Box& operator=(const Box<SourceScalarParam,sourceDimensionParam>& source) // Assignment with type conversion and dimension change
		{
		/* Operation is idempotent; no need to check for aliasing */
		min=source.min;
		max=source.max;
		return *this;
		}
	
	/* Methods: */
	friend bool operator==<>(const Box& b1,const Box& b2); // Equality operator
	friend bool operator!=<>(const Box& b1,const Box& b2); // Inequality operator
	bool isNull(void) const; // Returns true if the box contains no points at all
	bool isEmpty(void) const; // Returns true if the box has no interior (it can still contain points)
	bool isFull(void) const; // Returns true if the box contains all points
	const Point& getOrigin(void) const // Returns the box's origin
		{
		return min;
		}
	Box& setOrigin(const Point& newOrigin) // Sets a new origin while keeping the box's size
		{
		max=newOrigin+(max-min);
		min=newOrigin;
		return *this;
		}
	Size getSize(void) const // Returns the box's size
		{
		Size result;
		for(int i=0;i<dimension;++i)
			result[i]=max[i]-min[i];
		return result;
		}
	Box& setSize(const Size& newSize) // Sets a new size
		{
		for(int i=0;i<dimension;++i)
			max[i]=min[i]+newSize[i];
		return *this;
		}
	Scalar getSize(int index) const // Returns box's size along one primary axis
		{
		return max[index]-min[index];
		}
	Box& setSize(int index,Scalar newSize) // Sets a new size along one primary axis
		{
		max[index]=min[index]+newSize;
		return *this;
		}
	Point getVertex(int index) const // Returns one of the box's corner vertices (standard order)
		{
		Point result;
		for(int i=0;i<dimension;++i)
			result[i]=index&(1<<i)?max[i]:min[i];
		return result;
		}
	Box& setVertex(int index,const Point& newVertex) // Sets one of the box's corner vertices (standard order) to a new position
		{
		for(int i=0;i<dimension;++i)
			{
			if(index&(1<<i))
				max[i]=newVertex[i];
			else
				min[i]=newVertex[i];
			}
		return *this;
		}
	bool contains(const Point& p) const // Checks if a box contains a point
		{
		bool result=true;
		for(int i=0;i<dimension&&result;++i)
			result=min[i]<=p[i]&&p[i]<=max[i];
		return result;
		}
	bool contains(const Box& other) const // Checks if a box contains another box
		{
		bool result=true;
		for(int i=0;i<dimension&&result;++i)
			result=min[i]<=other.min[i]&&other.max[i]<=max[i];
		return result;
		}
	bool intersects(const Box& other) const // Checks if a box intersects another box
		{
		bool result=true;
		for(int i=0;i<dimension&&result;++i)
			result=min[i]<=other.max[i]&&other.min[i]<=max[i];
		return result;
		}
	bool overlaps(const Box& other) const // Checks if a box intersects another box with a non-zero volume intersection
		{
		bool result=true;
		for(int i=0;i<dimension&&result;++i)
			result=min[i]<other.max[i]&&other.min[i]<max[i];
		return result;
		}
	Box& shift(const Vector& offset) // Shifts a box by the given offset vector
		{
		min+=offset;
		max+=offset;
		return *this;
		}
	Box& extrude(Scalar d) // Extrudes all sides of a box by a given distance
		{
		for(int i=0;i<dimension;++i)
			{
			min[i]-=d;
			max[i]+=d;
			}
		return *this;
		}
	Box& extrude(const Size& s) // Extrudes box by given size in both directions
		{
		for(int i=0;i<dimension;++i)
			{
			min[i]-=s[i];
			max[i]+=s[i];
			}
		return *this;
		}
	Box& addPoint(const Point& p); // Changes the box to contain the given point
	Box& addBox(const Box& other); // Changes the box to contain the given box
	Box& intersectBox(const Box& other); // Changes the box to the intersection with the given box
	friend Box add<>(const Box& b1,const Box& b2); // Adds two boxes (union operator)
	friend Box intersect<>(const Box& b1,const Box& b2); // Intersects two boxes
	template <class TransformationParam>
	Box& transform(const TransformationParam& t) // Creates box containing the transformed interior of box
		{
		if(!(isNull()||isFull())) // Need not transform a box with no points inside; must not transform a box with all points inside
			{
			Point newMin,newMax;
			newMin=newMax=t.transform(min);
			for(int i=1;i<8;++i)
				{
				Point p=t.transform(getVertex(i));
				for(int j=0;j<dimension;++j)
					{
					if(newMin[j]>p[j])
						newMin[j]=p[j];
					else if(newMax[j]<p[j])
						newMax[j]=p[j];
					}
				}
			min=newMin;
			max=newMax;
			}
		return *this;
		}
	std::pair<Scalar,Scalar> getRayParameters(const Ray& ray) const; // Returns closed interval of ray parameters inside box
	HitResult intersectRay(const Ray& ray) const; // Intersects box with ray
	};

/*****************************
Friend functions of class Box:
*****************************/

template <class ScalarParam,int dimensionParam>
bool operator==(const Box<ScalarParam,dimensionParam>& b1,const Box<ScalarParam,dimensionParam>& b2)
	{
	return b1.min==b2.min&&b1.max==b2.max;
	}
template <class ScalarParam,int dimensionParam>
bool operator!=(const Box<ScalarParam,dimensionParam>& b1,const Box<ScalarParam,dimensionParam>& b2)
	{
	return b1.min!=b2.min||b1.max!=b2.max;
	}

}

#if defined(GEOMETRY_NONSTANDARD_TEMPLATES) && !defined(GEOMETRY_BOX_IMPLEMENTATION)
#include <Geometry/Box.icpp>
#endif

#endif
