/***********************************************************************
Vector - Class for Euclidean and affine vectors.
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

#ifndef GEOMETRY_VECTOR_INCLUDED
#define GEOMETRY_VECTOR_INCLUDED

#include <Math/Math.h>
#include <Geometry/ComponentArray.h>

namespace Geometry {

template <class ScalarParam,int dimensionParam>
class Vector:public ComponentArray<ScalarParam,dimensionParam>
	{
	/* Declarations of inherited types/elements: */
	public:
	using ComponentArray<ScalarParam,dimensionParam>::dimension;
	using ComponentArray<ScalarParam,dimensionParam>::components;
	
	/* Constructors and destructors: */
	public:
	static const Vector zero; // The zero vector
	Vector(void) // No initialization
		{
		}
	explicit Vector(ScalarParam filler) // Fills the component array with a single value
		:ComponentArray<ScalarParam,dimensionParam>(filler)
		{
		}
	Vector(ScalarParam sX,ScalarParam sY) // Constructor for 2D vector
		:ComponentArray<ScalarParam,dimensionParam>(sX,sY)
		{
		}
	Vector(ScalarParam sX,ScalarParam sY,ScalarParam sZ) // Constructor for 3D vector
		:ComponentArray<ScalarParam,dimensionParam>(sX,sY,sZ)
		{
		}
	Vector(ScalarParam sX,ScalarParam sY,ScalarParam sZ,ScalarParam sW) // Constructor for 4D vector
		:ComponentArray<ScalarParam,dimensionParam>(sX,sY,sZ,sW)
		{
		}
	template <class SourceScalarParam>
	Vector(const SourceScalarParam array[dimensionParam]) // Construction from C-style array
		:ComponentArray<ScalarParam,dimensionParam>(array)
		{
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	explicit Vector(const ComponentArray<SourceScalarParam,sourceDimensionParam>& source) // Constructs a vector from a component array
		:ComponentArray<ScalarParam,dimensionParam>(source)
		{
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	Vector(const Vector<SourceScalarParam,sourceDimensionParam>& source) // Copy constructor with type conversion and dimension change
		:ComponentArray<ScalarParam,dimensionParam>(source)
		{
		}
	
	/* Methods: */
	Vector operator+(void) const // Unary plus operator; returns copy of vector
		{
		return *this;
		}
	Vector operator-(void) const // Negation operator
		{
		Vector result;
		for(int i=0;i<dimension;++i)
			result.components[i]=-components[i];
		return result;
		}
	Vector& operator+=(const Vector& other) // Addition assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]+=other.components[i];
		return *this;
		}
	Vector& operator-=(const Vector& other) // Subtraction assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]-=other.components[i];
		return *this;
		}
	Vector& operator*=(ScalarParam scalar) // Scalar multiplication assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]*=scalar;
		return *this;
		}
	Vector& operator/=(ScalarParam scalar) // Scalar division assignment
		{
		for(int i=0;i<dimension;++i)
			components[i]/=scalar;
		return *this;
		}
	Vector& normalize(void) // Scales a vector to unit length
		{
		double norm=0.0;
		for(int i=0;i<dimension;++i)
			norm+=Math::sqr(double(components[i]));
		norm=Math::sqrt(norm);
		for(int i=0;i<dimension;++i)
			components[i]/=norm;
		return *this;
		}
	Vector& orthogonalize(const Vector& normal) // Orthogonalizes a vector with respect to the (non-unit length) normal vector
		{
		ScalarParam proj(0);
		ScalarParam denom(0);
		for(int i=0;i<dimension;++i)
			{
			proj+=components[i]*normal.components[i];
			denom+=normal.components[i]*normal.components[i];
			}
		proj/=denom;
		for(int i=0;i<dimension;++i)
			components[i]-=normal.components[i]*proj;
		return *this;
		}
	Vector& reflect(const Vector& normal) // Reflects a vector with respect to the plane defined by the (non-unit length) normal vector
		{
		ScalarParam proj(0);
		ScalarParam denom(0);
		for(int i=0;i<dimension;++i)
			{
			proj+=components[i]*normal.components[i];
			denom+=normal.components[i]*normal.components[i];
			}
		proj=ScalarParam(2)*proj/denom;
		for(int i=0;i<dimension;++i)
			components[i]=normal.components[i]*proj-components[i];
		return *this;
		}
	// Vector cross(const Vector& other) const; // Returns cross product of two vectors
	};

#if 0
template<>
inline Vector<int,3> Vector<int,3>::cross(const Vector<int,3>& other) const
	{
	return Vector(components[1]*other.components[2]-components[2]*other.components[1],components[2]*other.components[0]-components[0]*other.components[2],components[0]*other.components[1]-components[1]*other.components[0]);
	}

template<>
inline Vector<float,3> Vector<float,3>::cross(const Vector<float,3>& other) const
	{
	return Vector(components[1]*other.components[2]-components[2]*other.components[1],components[2]*other.components[0]-components[0]*other.components[2],components[0]*other.components[1]-components[1]*other.components[0]);
	}

template<>
inline Vector<double,3> Vector<double,3>::cross(const Vector<double,3>& other) const
	{
	return Vector(components[1]*other.components[2]-components[2]*other.components[1],components[2]*other.components[0]-components[0]*other.components[2],components[0]*other.components[1]-components[1]*other.components[0]);
	}
#endif

/*************************************
Operations on objects of class Vector:
*************************************/

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator+(const Vector<ScalarParam,dimensionParam>& v1,const Vector<ScalarParam,dimensionParam>& v2) // Addition
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=v1[i]+v2[i];
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator+(const Vector<ScalarParam,2>& v1,const Vector<ScalarParam,2>& v2)
	{
	return Vector<ScalarParam,2>(v1[0]+v2[0],v1[1]+v2[1]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator+(const Vector<ScalarParam,3>& v1,const Vector<ScalarParam,3>& v2)
	{
	return Vector<ScalarParam,3>(v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,4> operator+(const Vector<ScalarParam,4>& v1,const Vector<ScalarParam,4>& v2)
	{
	return Vector<ScalarParam,4>(v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2],v1[3]+v2[3]);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator-(const Vector<ScalarParam,dimensionParam>& v1,const Vector<ScalarParam,dimensionParam>& v2) // Subtraction
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=v1[i]-v2[i];
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator-(const Vector<ScalarParam,2>& v1,const Vector<ScalarParam,2>& v2)
	{
	return Vector<ScalarParam,2>(v1[0]-v2[0],v1[1]-v2[1]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator-(const Vector<ScalarParam,3>& v1,const Vector<ScalarParam,3>& v2)
	{
	return Vector<ScalarParam,3>(v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,4> operator-(const Vector<ScalarParam,4>& v1,const Vector<ScalarParam,4>& v2)
	{
	return Vector<ScalarParam,4>(v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2],v1[3]-v2[3]);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator*(const Vector<ScalarParam,dimensionParam>& v,ScalarParam scalar) // Scalar multiplication (from the right)
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=v[i]*scalar;
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator*(const Vector<ScalarParam,2>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,2>(v[0]*scalar,v[1]*scalar);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator*(const Vector<ScalarParam,3>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,3>(v[0]*scalar,v[1]*scalar,v[2]*scalar);
	}

template <class ScalarParam>
inline Vector<ScalarParam,4> operator*(const Vector<ScalarParam,4>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,4>(v[0]*scalar,v[1]*scalar,v[2]*scalar,v[3]*scalar);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator*(ScalarParam scalar,const Vector<ScalarParam,dimensionParam>& v1) // Ditto (from the left)
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=scalar*v1[i];
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator*(ScalarParam scalar,const Vector<ScalarParam,2>& v)
	{
	return Vector<ScalarParam,2>(scalar*v[0],scalar*v[1]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator*(ScalarParam scalar,const Vector<ScalarParam,3>& v)
	{
	return Vector<ScalarParam,3>(scalar*v[0],scalar*v[1],scalar*v[2]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,4> operator*(ScalarParam scalar,const Vector<ScalarParam,4>& v)
	{
	return Vector<ScalarParam,4>(scalar*v[0],scalar*v[1],scalar*v[2],scalar*v[3]);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> operator/(const Vector<ScalarParam,dimensionParam>& v1,ScalarParam scalar) // Scalar division
	{
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=v1[i]/scalar;
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> operator/(const Vector<ScalarParam,2>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,2>(v[0]/scalar,v[1]/scalar);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator/(const Vector<ScalarParam,3>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,3>(v[0]/scalar,v[1]/scalar,v[2]/scalar);
	}

template <class ScalarParam>
inline Vector<ScalarParam,4> operator/(const Vector<ScalarParam,4>& v,ScalarParam scalar)
	{
	return Vector<ScalarParam,4>(v[0]/scalar,v[1]/scalar,v[2]/scalar,v[3]/scalar);
	}

template <class ScalarParam,int dimensionParam>
Vector<ScalarParam,dimensionParam> normalize(const Vector<ScalarParam,dimensionParam>& v) // Returns a collinear vector of unit length
	{
	double norm=0.0;
	for(int i=0;i<dimensionParam;++i)
		norm+=Math::sqr(double(v[i]));
	norm=Math::sqrt(norm);
	Vector<ScalarParam,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[i]=ScalarParam(v[i]/norm);
	return result;
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> orthogonalize(const Vector<ScalarParam,dimensionParam>& v,const Vector<ScalarParam,dimensionParam>& normal) // Orthogonalizes vector v with respect to the (non-unit length) normal vector
	{
	return v-normal*((v*normal)/normal.sqr());
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> cross(const Vector<ScalarParam,3>& v1,const Vector<ScalarParam,3>& v2) // Returns the cross product of two 3D vectors
	{
	return Vector<ScalarParam,3>(v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]);
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> operator^(const Vector<ScalarParam,3>& v1,const Vector<ScalarParam,3>& v2) // Returns the cross product of two 3D vectors
	{
	return Vector<ScalarParam,3>(v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]);
	}

template <class ScalarParam,int dimensionParam>
inline Vector<ScalarParam,dimensionParam> reflect(const Vector<ScalarParam,dimensionParam>& v,const Vector<ScalarParam,dimensionParam>& normal) // Reflects vector v with respect to the plane orthogonal to the (non-unit length) normal vector
	{
	return normal*(ScalarParam(2)*(v*normal)/normal.sqr())-v;
	}

template <class ScalarParam,int dimensionParam>
inline int findParallelAxis(const Vector<ScalarParam,dimensionParam>& v) // Finds the index of the primary axis most parallel to the vector
	{
	int result=0;
	ScalarParam maxAbsAxis=Math::abs(v[0]);
	for(int i=1;i<dimensionParam;++i)
		{
		ScalarParam absAxis=Math::abs(v[i]);
		if(maxAbsAxis<absAxis)
			{
			result=i;
			maxAbsAxis=absAxis;
			}
		}
	return result;
	}

template <class ScalarParam,int dimensionParam>
inline int findOrthogonalAxis(const Vector<ScalarParam,dimensionParam>& v) // Finds the index of the primary axis most orthogonal to the vector
	{
	int result=0;
	ScalarParam minAbsAxis=Math::abs(v[0]);
	for(int i=1;i<dimensionParam;++i)
		{
		ScalarParam absAxis=Math::abs(v[i]);
		if(minAbsAxis>absAxis)
			{
			result=i;
			minAbsAxis=absAxis;
			}
		}
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,2> normal(const Vector<ScalarParam,2>& v)
	{
	Vector<ScalarParam,2> result;
	result[0]=-v[1];
	result[1]=v[0];
	return result;
	}

template <class ScalarParam>
inline Vector<ScalarParam,3> normal(const Vector<ScalarParam,3>& v)
	{
	ScalarParam t[3];
	for(int i=0;i<3;++i)
		t[i]=Math::abs(v[i]);
	Vector<ScalarParam,3> result;
	if(t[0]<t[1]&&t[0]<t[2])
		{
		result[0]=ScalarParam(0);
		result[1]=v[2];
		result[2]=-v[1];
		}
	else if(t[1]<t[2])
		{
		result[0]=v[2];
		result[1]=ScalarParam(0);
		result[2]=-v[0];
		}
	else
		{
		result[0]=v[1];
		result[1]=-v[0];
		result[2]=ScalarParam(0);
		}
	return result;
	}

}

#if defined(GEOMETRY_NONSTANDARD_TEMPLATES) && !defined(GEOMETRY_VECTOR_IMPLEMENTATION)
#include <Geometry/Vector.icpp>
#endif

#endif
