/***********************************************************************
ComponentArray - Foundation class for Euclidean vectors, affine points
and homogenuous vectors.
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

#ifndef GEOMETRY_COMPONENTARRAY_INCLUDED
#define GEOMETRY_COMPONENTARRAY_INCLUDED

#include <Math/Math.h>

namespace Geometry {

/********************************
The generic ComponentArray class:
********************************/

template <class ScalarParam,int dimensionParam>
class ComponentArray
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=dimensionParam; // The component array's dimension
	
	/* Elements: */
	protected:
	Scalar components[dimension]; // The component array
	
	/* Constructors and destructors: */
	public:
	ComponentArray(void) // No initialization
		{
		}
	explicit ComponentArray(Scalar filler) // Fills the component array with a single value
		{
		for(int i=0;i<dimension;++i)
			components[i]=filler;
		}
	template <class SourceScalarParam>
	ComponentArray(const SourceScalarParam array[dimensionParam]) // Construction from C-style array
		{
		/* Copy and typecast the source array: */
		for(int i=0;i<dimension;++i)
			components[i]=Scalar(array[i]);
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	ComponentArray(const ComponentArray<SourceScalarParam,sourceDimensionParam>& source); // Copy constructor with type conversion and dimension change
	template <class SourceScalarParam>
	ComponentArray(const ComponentArray<SourceScalarParam,dimensionParam>& source); // Copy constructor with type conversion
	
	/* Methods: */
	const Scalar* getComponents(void) const // Returns component array
		{
		return components;
		}
	Scalar* getComponents(void) // Ditto
		{
		return components;
		}
	Scalar operator[](int index) const // Returns component as rvalue
		{
		return components[index];
		}
	Scalar& operator[](int index) // Returns component as modifiable lvalue
		{
		return components[index];
		}
	
	/* Some common operations on ComponentArrays: */
	Scalar sqr(void) const // Returns the squared L2 norm of a component array
		{
		Scalar result(0);
		for(int i=0;i<dimension;++i)
			result+=Math::sqr(components[i]);
		return result;
		}
	double mag(void) const // Returns the L2 norm of a component array
		{
		double result=0.0;
		for(int i=0;i<dimension;++i)
			result+=Math::sqr(double(components[i]));
		return Math::sqrt(result);
		}
	Scalar abs(void) const // Returns the absolute norm of a component array
		{
		Scalar result(0);
		for(int i=0;i<dimension;++i)
			result+=Math::abs(components[i]);
		return result;
		}
	Scalar max(void) const // Returns the maximum norm of a component array
		{
		Scalar result=Math::abs(components[0]);
		for(int i=1;i<dimension;++i)
			{
			Scalar val=Math::abs(components[i]);
			if(result<val)
				result=val;
			}
		return result;
		}
	};

/************************************************
Specialized versions of the ComponentArray class:
************************************************/

template <class ScalarParam>
class ComponentArray<ScalarParam,2>
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=2; // The component array's dimension
	
	/* Elements: */
	protected:
	Scalar components[dimension]; // The component array
	
	/* Constructors and destructors: */
	public:
	ComponentArray(void) // No initialization
		{
		}
	explicit ComponentArray(Scalar filler) // Fills the component array with a single value
		{
		components[0]=filler;
		components[1]=filler;
		}
	ComponentArray(Scalar c0,Scalar c1)
		{
		components[0]=c0;
		components[1]=c1;
		}
	template <class SourceScalarParam>
	ComponentArray(const SourceScalarParam array[2]) // Construction from C-style array
		{
		components[0]=Scalar(array[0]);
		components[1]=Scalar(array[1]);
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	ComponentArray(const ComponentArray<SourceScalarParam,sourceDimensionParam>& source); // Copy constructor with type conversion and dimension change
	template <class SourceScalarParam>
	ComponentArray(const ComponentArray<SourceScalarParam,2>& source); // Copy constructor with type conversion
	
	/* Methods: */
	const Scalar* getComponents(void) const // Returns component array
		{
		return components;
		}
	Scalar* getComponents(void) // Ditto
		{
		return components;
		}
	Scalar operator[](int index) const // Returns component as rvalue
		{
		return components[index];
		}
	Scalar& operator[](int index) // Returns component as modifiable lvalue
		{
		return components[index];
		}
	
	/* Some common operations on ComponentArrays: */
	Scalar sqr(void) const // Returns the squared L2 norm of a component array
		{
		return Math::sqr(components[0])+Math::sqr(components[1]);
		}
	double mag(void) const // Returns the L2 norm of a component array
		{
		return Math::sqrt(Math::sqr(double(components[0]))+Math::sqr(double(components[1])));
		}
	Scalar abs(void) const // Returns the absolute norm of a component array
		{
		return Math::abs(components[0])+Math::abs(components[1]);
		}
	Scalar max(void) const // Returns the maximum norm of a component array
		{
		Scalar result=Math::abs(components[0]);
		Scalar val;
		if(result<(val=Math::abs(components[1])))
			result=val;
		return result;
		}
	};

template <class ScalarParam>
class ComponentArray<ScalarParam,3>
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=3; // The component array's dimension
	
	/* Elements: */
	protected:
	Scalar components[dimension]; // The component array
	
	/* Constructors and destructors: */
	public:
	ComponentArray(void) // No initialization
		{
		}
	explicit ComponentArray(Scalar filler) // Fills the component array with a single value
		{
		components[0]=filler;
		components[1]=filler;
		components[2]=filler;
		}
	ComponentArray(Scalar c0,Scalar c1,Scalar c2 =Scalar(0))
		{
		components[0]=c0;
		components[1]=c1;
		components[2]=c2;
		}
	template <class SourceScalarParam>
	ComponentArray(const SourceScalarParam array[3]) // Construction from C-style array
		{
		components[0]=Scalar(array[0]);
		components[1]=Scalar(array[1]);
		components[2]=Scalar(array[2]);
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	ComponentArray(const ComponentArray<SourceScalarParam,sourceDimensionParam>& source); // Copy constructor with type conversion and dimension change
	template <class SourceScalarParam>
	ComponentArray(const ComponentArray<SourceScalarParam,3>& source); // Copy constructor with type conversion
	
	/* Methods: */
	const ScalarParam* getComponents(void) const // Returns component array
		{
		return components;
		}
	Scalar* getComponents(void) // Ditto
		{
		return components;
		}
	Scalar operator[](int index) const // Returns component as rvalue
		{
		return components[index];
		}
	Scalar& operator[](int index) // Returns component as modifiable lvalue
		{
		return components[index];
		}
	
	/* Some common operations on ComponentArrays: */
	Scalar sqr(void) const // Returns the squared L2 norm of a component array
		{
		return Math::sqr(components[0])+Math::sqr(components[1])+Math::sqr(components[2]);
		}
	double mag(void) const // Returns the L2 norm of a component array
		{
		return Math::sqrt(Math::sqr(double(components[0]))+Math::sqr(double(components[1]))+Math::sqr(double(components[2])));
		}
	Scalar abs(void) const // Returns the absolute norm of a component array
		{
		return Math::abs(components[0])+Math::abs(components[1])+Math::abs(components[2]);
		}
	Scalar max(void) const // Returns the maximum norm of a component array
		{
		Scalar result=Math::abs(components[0]);
		Scalar val;
		if(result<(val=Math::abs(components[1])))
			result=val;
		if(result<(val=Math::abs(components[2])))
			result=val;
		return result;
		}
	};

template <class ScalarParam>
class ComponentArray<ScalarParam,4>
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar; // The underlying scalar type
	static const int dimension=4; // The component array's dimension
	
	/* Elements: */
	protected:
	Scalar components[dimension]; // The component array
	
	/* Constructors and destructors: */
	public:
	ComponentArray(void) // No initialization
		{
		}
	explicit ComponentArray(Scalar filler) // Fills the component array with a single value
		{
		for(int i=0;i<4;++i)
			components[i]=filler;
		}
	ComponentArray(Scalar c0,Scalar c1,Scalar c2 =Scalar(0),Scalar c3 =Scalar(0))
		{
		components[0]=c0;
		components[1]=c1;
		components[2]=c2;
		components[3]=c3;
		}
	template <class SourceScalarParam>
	ComponentArray(const SourceScalarParam array[4]) // Construction from C-style array
		{
		for(int i=0;i<4;++i)
			components[i]=Scalar(array[i]);
		}
	template <class SourceScalarParam,int sourceDimensionParam>
	ComponentArray(const ComponentArray<SourceScalarParam,sourceDimensionParam>& source); // Copy constructor with type conversion and dimension change
	template <class SourceScalarParam>
	ComponentArray(const ComponentArray<SourceScalarParam,4>& source); // Copy constructor with type conversion
	
	/* Methods: */
	const Scalar* getComponents(void) const // Returns component array
		{
		return components;
		}
	Scalar* getComponents(void) // Ditto
		{
		return components;
		}
	Scalar operator[](int index) const // Returns component as rvalue
		{
		return components[index];
		}
	Scalar& operator[](int index) // Returns component as modifiable lvalue
		{
		return components[index];
		}
	
	/* Some common operations on ComponentArrays: */
	Scalar sqr(void) const // Returns the squared L2 norm of a component array
		{
		return Math::sqr(components[0])+Math::sqr(components[1])+Math::sqr(components[2])+Math::sqr(components[3]);
		}
	double mag(void) const // Returns the L2 norm of a component array
		{
		return Math::sqrt(Math::sqr(double(components[0]))+Math::sqr(double(components[1]))+Math::sqr(double(components[2]))+Math::sqr(double(components[3])));
		}
	Scalar abs(void) const // Returns the absolute norm of a component array
		{
		return Math::abs(components[0])+Math::abs(components[1])+Math::abs(components[2])+Math::abs(components[3]);
		}
	Scalar max(void) const // Returns the maximum norm of a component array
		{
		Scalar result=Math::abs(components[0]);
		for(int i=1;i<dimension;++i)
			{
			Scalar val=Math::abs(components[i]);
			if(result<val)
				result=val;
			}
		return result;
		}
	};

/*********************************************
Operations on objects of class ComponentArray:
*********************************************/

template <class ScalarParam,int dimensionParam>
inline bool operator==(const ComponentArray<ScalarParam,dimensionParam>& ca1,const ComponentArray<ScalarParam,dimensionParam>& ca2) // Equality operator
	{
	bool result=true;
	for(int i=0;i<dimensionParam;++i)
		result&=ca1[i]==ca2[i];
	return result;
	}

template <class ScalarParam>
inline bool operator==(const ComponentArray<ScalarParam,2>& ca1,const ComponentArray<ScalarParam,2>& ca2)
	{
	return ca1[0]==ca2[0]&&ca1[1]==ca2[1];
	}

template <class ScalarParam>
inline bool operator==(const ComponentArray<ScalarParam,3>& ca1,const ComponentArray<ScalarParam,3>& ca2)
	{
	return ca1[0]==ca2[0]&&ca1[1]==ca2[1]&&ca1[2]==ca2[2];
	}

template <class ScalarParam>
inline bool operator==(const ComponentArray<ScalarParam,4>& ca1,const ComponentArray<ScalarParam,4>& ca2)
	{
	return ca1[0]==ca2[0]&&ca1[1]==ca2[1]&&ca1[2]==ca2[2]&&ca1[3]==ca2[3];
	}

template <class ScalarParam,int dimensionParam>
inline bool operator!=(const ComponentArray<ScalarParam,dimensionParam>& ca1,const ComponentArray<ScalarParam,dimensionParam>& ca2) // Inequality operator
	{
	bool result=false;
	for(int i=0;i<dimensionParam;++i)
		result|=ca1[i]!=ca2[i];
	return result;
	}

template <class ScalarParam>
inline bool operator!=(const ComponentArray<ScalarParam,2>& ca1,const ComponentArray<ScalarParam,2>& ca2)
	{
	return ca1[0]!=ca2[0]||ca1[1]!=ca2[1];
	}

template <class ScalarParam>
inline bool operator!=(const ComponentArray<ScalarParam,3>& ca1,const ComponentArray<ScalarParam,3>& ca2)
	{
	return ca1[0]!=ca2[0]||ca1[1]!=ca2[1]||ca1[2]!=ca2[2];
	}

template <class ScalarParam>
inline bool operator!=(const ComponentArray<ScalarParam,4>& ca1,const ComponentArray<ScalarParam,4>& ca2)
	{
	return ca1[0]!=ca2[0]||ca1[1]!=ca2[1]||ca1[2]!=ca2[2]||ca1[3]!=ca2[3];
	}

template <class ScalarParam,int dimensionParam>
inline ScalarParam operator*(const ComponentArray<ScalarParam,dimensionParam>& ca1,const ComponentArray<ScalarParam,dimensionParam>& ca2) // Returns the inner product of two component arrays
	{
	ScalarParam result(0);
	for(int i=0;i<dimensionParam;++i)
		result+=ca1[i]*ca2[i];
	return result;
	}

template <class ScalarParam>
inline ScalarParam operator*(const ComponentArray<ScalarParam,2>& ca1,const ComponentArray<ScalarParam,2>& ca2) // Ditto
	{
	return ca1[0]*ca2[0]+ca1[1]*ca2[1];
	}

template <class ScalarParam>
inline ScalarParam operator*(const ComponentArray<ScalarParam,3>& ca1,const ComponentArray<ScalarParam,3>& ca2) // Ditto
	{
	return ca1[0]*ca2[0]+ca1[1]*ca2[1]+ca1[2]*ca2[2];
	}

template <class ScalarParam>
inline ScalarParam operator*(const ComponentArray<ScalarParam,4>& ca1,const ComponentArray<ScalarParam,4>& ca2) // Ditto
	{
	return ca1[0]*ca2[0]+ca1[1]*ca2[1]+ca1[2]*ca2[2]+ca1[3]*ca2[3];
	}

template <class ScalarParam,int dimensionParam>
inline ScalarParam sqr(const ComponentArray<ScalarParam,dimensionParam>& ca) // Returns the squared L2 norm of a component array
	{
	ScalarParam result(0);
	for(int i=0;i<dimensionParam;++i)
		result+=Math::sqr(ca[i]);
	return result;
	}

template <class ScalarParam>
inline ScalarParam sqr(const ComponentArray<ScalarParam,2>& ca) // Ditto
	{
	return Math::sqr(ca[0])+Math::sqr(ca[1]);
	}

template <class ScalarParam>
inline ScalarParam sqr(const ComponentArray<ScalarParam,3>& ca) // Ditto
	{
	return Math::sqr(ca[0])+Math::sqr(ca[1])+Math::sqr(ca[2]);
	}

template <class ScalarParam>
inline ScalarParam sqr(const ComponentArray<ScalarParam,4>& ca) // Ditto
	{
	return Math::sqr(ca[0])+Math::sqr(ca[1])+Math::sqr(ca[2])+Math::sqr(ca[3]);
	}

template <class ScalarParam,int dimensionParam>
inline double mag(const ComponentArray<ScalarParam,dimensionParam>& ca) // Returns the L2 norm of a component array
	{
	double result=0.0;
	for(int i=0;i<dimensionParam;++i)
		result+=Math::sqr(double(ca[i]));
	return Math::sqrt(result);
	}

template <class ScalarParam,int dimensionParam>
inline ScalarParam abs(const ComponentArray<ScalarParam,dimensionParam>& ca) // Returns the absolute norm of a component array
	{
	ScalarParam result(0);
	for(int i=0;i<dimensionParam;++i)
		result+=Math::abs(ca[i]);
	return result;
	}

template <class ScalarParam,int dimensionParam>
inline ScalarParam max(const ComponentArray<ScalarParam,dimensionParam>& ca) // Returns the maximum norm of a component array
	{
	ScalarParam result=Math::abs(ca[0]);
	for(int i=1;i<dimensionParam;++i)
		{
		ScalarParam val=Math::abs(ca[i]);
		if(result<val)
			result=val;
		}
	return result;
	}

}

#if defined(GEOMETRY_NONSTANDARD_TEMPLATES) && !defined(GEOMETRY_COMPONENTARRAY_IMPLEMENTATION)
#include <Geometry/ComponentArray.icpp>
#endif

#endif
