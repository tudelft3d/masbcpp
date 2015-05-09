/***********************************************************************
Matrix - Class for n x m matrices, used internally in the implementation
of affine and projective transformations.
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

#ifndef GEOMETRY_MATRIX_INCLUDED
#define GEOMETRY_MATRIX_INCLUDED

#include <Geometry/ComponentArray.h>
#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/HVector.h>

namespace Geometry {

/* Helper class to specialize matrix operations: */

template <class ScalarParam,int numRowsParam,int numColumnsParam>
class MatrixOperations
	{
	/* Methods: */
	public:
	inline static ComponentArray<ScalarParam,numRowsParam> multiply(const ScalarParam m[numRowsParam][numColumnsParam],const ComponentArray<ScalarParam,numColumnsParam>& ca)
		{
		ComponentArray<ScalarParam,numRowsParam> result;
		for(int i=0;i<numRowsParam;++i)
			{
			ScalarParam comp(0);
			for(int j=0;j<numColumnsParam;++j)
				comp+=m[i][j]*ca[j];
			result[i]=comp;
			}
		return result;
		}
	inline static ComponentArray<ScalarParam,numColumnsParam> transposeMultiply(const ScalarParam m[numRowsParam][numColumnsParam],const ComponentArray<ScalarParam,numRowsParam>& ca)
		{
		ComponentArray<ScalarParam,numColumnsParam> result;
		for(int j=0;j<numColumnsParam;++j)
			{
			ScalarParam comp(0);
			for(int i=0;i<numRowsParam;++i)
				comp+=m[i][j]*ca[i];
			result[j]=comp;
			}
		return result;
		}
	};

template <class ScalarParam,int dimensionParam>
class MatrixOperations<ScalarParam,dimensionParam,dimensionParam>
	{
	/* Methods: */
	public:
	inline static ComponentArray<ScalarParam,dimensionParam> multiply(const ScalarParam m[dimensionParam][dimensionParam],const ComponentArray<ScalarParam,dimensionParam>& ca)
		{
		ComponentArray<ScalarParam,dimensionParam> result;
		for(int i=0;i<dimensionParam;++i)
			{
			ScalarParam comp(0);
			for(int j=0;j<dimensionParam;++j)
				comp+=m[i][j]*ca[j];
			result[i]=comp;
			}
		return result;
		}
	inline static ComponentArray<ScalarParam,dimensionParam> transposeMultiply(const ScalarParam m[dimensionParam][dimensionParam],const ComponentArray<ScalarParam,dimensionParam>& ca)
		{
		ComponentArray<ScalarParam,dimensionParam> result;
		for(int j=0;j<dimensionParam;++j)
			{
			ScalarParam comp(0);
			for(int i=0;i<dimensionParam;++i)
				comp+=m[i][j]*ca[i];
			result[j]=comp;
			}
		return result;
		}
	static ComponentArray<ScalarParam,dimensionParam> divide(const ComponentArray<ScalarParam,dimensionParam>& ca,const ScalarParam m[dimensionParam][dimensionParam]);
	static double determinant(const ScalarParam m[dimensionParam][dimensionParam]);
	static void invert(const ScalarParam m[dimensionParam][dimensionParam],ScalarParam result[dimensionParam][dimensionParam]);
	};

template <class ScalarParam>
class MatrixOperations<ScalarParam,2,2>
	{
	/* Methods: */
	public:
	inline static ComponentArray<ScalarParam,2> multiply(const ScalarParam m[2][2],const ComponentArray<ScalarParam,2>& ca)
		{
		return ComponentArray<ScalarParam,2>(m[0][0]*ca[0]+m[0][1]*ca[1],
		                                     m[1][0]*ca[0]+m[1][1]*ca[1]);
		}
	inline static ComponentArray<ScalarParam,2> transposeMultiply(const ScalarParam m[2][2],const ComponentArray<ScalarParam,2>& ca)
		{
		return ComponentArray<ScalarParam,2>(m[0][0]*ca[0]+m[1][0]*ca[1],
		                                     m[0][1]*ca[0]+m[1][1]*ca[1]);
		}
	inline static ComponentArray<ScalarParam,2> divide(const ComponentArray<ScalarParam,2>& ca,const ScalarParam m[2][2])
		{
		double det=double(m[0][0])*double(m[1][1])-double(m[1][0])*double(m[0][1]);
		return ComponentArray<ScalarParam,2>(ScalarParam((double(m[1][1])*double(ca[0])-double(m[0][1])*double(ca[1]))/det),
		                                     ScalarParam((double(m[0][0])*double(ca[1])-double(m[1][0])*double(ca[0]))/det));
		}
	inline static double determinant(const ScalarParam m[2][2])
		{
		return double(m[0][0])*double(m[1][1])-double(m[1][0])*double(m[0][1]);
		}
	inline static void invert(const ScalarParam m[2][2],ScalarParam result[2][2])
		{
		double det=double(m[0][0])*double(m[1][1])-double(m[1][0])*double(m[0][1]);
		ScalarParam r11=ScalarParam(m[0][0]/det);
		result[0][0]=ScalarParam(m[1][1]/det);
		result[0][1]=ScalarParam(-m[0][1]/det);
		result[1][0]=ScalarParam(-m[1][0]/det);
		result[1][1]=r11;
		}
	};

template <class ScalarParam>
class MatrixOperations<ScalarParam,3,3>
	{
	/* Methods: */
	public:
	inline static ComponentArray<ScalarParam,3> multiply(const ScalarParam m[3][3],const ComponentArray<ScalarParam,3>& ca)
		{
		return ComponentArray<ScalarParam,3>(m[0][0]*ca[0]+m[0][1]*ca[1]+m[0][2]*ca[2],
		                                     m[1][0]*ca[0]+m[1][1]*ca[1]+m[1][2]*ca[2],
		                                     m[2][0]*ca[0]+m[2][1]*ca[1]+m[2][2]*ca[2]);
		}
	inline static ComponentArray<ScalarParam,3> transposeMultiply(const ScalarParam m[3][3],const ComponentArray<ScalarParam,3>& ca)
		{
		return ComponentArray<ScalarParam,3>(m[0][0]*ca[0]+m[1][0]*ca[1]+m[2][0]*ca[2],
		                                     m[0][1]*ca[0]+m[1][1]*ca[1]+m[2][1]*ca[2],
		                                     m[0][2]*ca[0]+m[1][2]*ca[1]+m[2][2]*ca[2]);
		}
	inline static ComponentArray<ScalarParam,3> divide(const ComponentArray<ScalarParam,3>& ca,const ScalarParam m[3][3])
		{
		double sub[3][3];
		sub[0][0]=double(m[1][1])*double(m[2][2])-double(m[2][1])*double(m[1][2]);
		sub[0][1]=double(m[1][2])*double(m[2][0])-double(m[2][2])*double(m[1][0]);
		sub[0][2]=double(m[1][0])*double(m[2][1])-double(m[2][0])*double(m[1][1]);
		sub[1][0]=double(m[2][1])*double(m[0][2])-double(m[0][1])*double(m[2][2]);
		sub[1][1]=double(m[2][2])*double(m[0][0])-double(m[0][2])*double(m[2][0]);
		sub[1][2]=double(m[2][0])*double(m[0][1])-double(m[0][0])*double(m[2][1]);
		sub[2][0]=double(m[0][1])*double(m[1][2])-double(m[1][1])*double(m[0][2]);
		sub[2][1]=double(m[0][2])*double(m[1][0])-double(m[1][2])*double(m[0][0]);
		sub[2][2]=double(m[0][0])*double(m[1][1])-double(m[1][0])*double(m[0][1]);
		double det=double(m[0][0])*sub[0][0]+double(m[1][0])*sub[1][0]+double(m[2][0])*sub[2][0];
		return ComponentArray<ScalarParam,3>(ScalarParam((sub[0][0]*double(ca[0])+sub[1][0]*double(ca[1])+sub[2][0]*double(ca[2]))/det),
		                                     ScalarParam((sub[0][1]*double(ca[0])+sub[1][1]*double(ca[1])+sub[2][1]*double(ca[2]))/det),
		                                     ScalarParam((sub[0][2]*double(ca[0])+sub[1][2]*double(ca[1])+sub[2][2]*double(ca[2]))/det));
		}
	inline static double determinant(const ScalarParam m[3][3])
		{
		return double(m[0][0])*(double(m[1][1])*double(m[2][2])-double(m[2][1])*double(m[1][2]))+
		       double(m[1][0])*(double(m[2][1])*double(m[0][2])-double(m[0][1])*double(m[2][2]))+
		       double(m[2][0])*(double(m[0][1])*double(m[1][2])-double(m[1][1])*double(m[0][2]));
		}
	static void invert(const ScalarParam m[3][3],ScalarParam result[3][3]);
	};

/* Forward declarations: */
template <class ScalarParam,int numRowsParam,int numColumnsParam>
class Matrix;

/* Forward declarations for friend functions: */
template <class ScalarParam,int numRowsParam,int numColumnsParam>
ComponentArray<ScalarParam,numRowsParam> operator*(const Matrix<ScalarParam,numRowsParam,numColumnsParam>&,const ComponentArray<ScalarParam,numColumnsParam>&);
template <class ScalarParam,int numRowsParam,int numColumnsParam>
ComponentArray<ScalarParam,numColumnsParam> operator/(const ComponentArray<ScalarParam,numRowsParam>&,const Matrix<ScalarParam,numRowsParam,numColumnsParam>&);
template <class ScalarParam,int numRowsParam,int numColumnsParam>
double determinant(const Matrix<ScalarParam,numRowsParam,numColumnsParam>&);
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numColumnsParam,numRowsParam> invert(const Matrix<ScalarParam,numRowsParam,numColumnsParam>&);

template <class ScalarParam,int numRowsParam,int numColumnsParam>
class Matrix
	{
	/* Embedded classes: */
	public:
	typedef ScalarParam Scalar;
	static const int numRows=numRowsParam;
	static const int numColumns=numColumnsParam;
	private:
	typedef MatrixOperations<ScalarParam,numRowsParam,numColumnsParam> MO;
	
	/* Elements: */
	private:
	Scalar c[numRows][numColumns]; // Array of matrix entries in row-major order
	
	/* Constructors and destructors: */
	public:
	static const Matrix zero; // The zero matrix
	static const Matrix one; // The unity matrix
	Matrix(void) // Creates uninitialized matrix
		{
		}
	Matrix(ScalarParam diagonal) // Creates matrix from diagonal element
		{
		for(int i=0;i<numRows;++i)
			for(int j=0;j<numColumns;++j)
				c[i][j]=i==j?diagonal:Scalar(0);
		}
	template <class SourceScalarParam>
	static Matrix fromRowMajor(const SourceScalarParam* components); // Constructs matrix from array in row-major order with type conversion
	template <class SourceScalarParam>
	static Matrix fromColumnMajor(const SourceScalarParam* components); // Constructs matrix from array in column-major order with type conversion
	template <class SourceScalarParam>
	Matrix(const Matrix<SourceScalarParam,numRowsParam,numColumnsParam>& source); // Copy constructor with type conversion
	
	/* Methods: */
	const Scalar* getEntries(void) const // Returns pointer to array of matrix entries
		{
		return &c[0][0];
		}
	Scalar* getEntries(void) // Ditto
		{
		return &c[0][0];
		}
	Scalar operator()(int row,int column) const // Returns a matrix entry as rvalue
		{
		return c[row][column];
		}
	Scalar& operator()(int row,int column) // Ditto, as modifiable lvalue
		{
		return c[row][column];
		}
	Matrix operator+(void) const // Unary plus returns copy of matrix
		{
		return *this;
		}
	Matrix operator-(void) const; // Negation operator
	Matrix& operator+=(const Matrix& other); // Addition assignment
	Matrix& operator-=(const Matrix& other); // Subtraction assignment
	Matrix& operator*=(Scalar scalar); // Scalar multiplication assignment
	Matrix& operator*=(const Matrix<Scalar,numColumnsParam,numColumnsParam>& other); // Matrix multiplication assignment (from the right)
	Matrix& transposeMultiply(const Matrix<Scalar,numColumnsParam,numColumnsParam>& other); // Matrix multiplication assignment with transposed matrix (from the right)
	Matrix& leftMultiply(const Matrix<Scalar,numRowsParam,numRowsParam>& other); // Matrix multiplication assignment (from the left)
	Matrix& transposeLeftMultiply(const Matrix<Scalar,numRowsParam,numRowsParam>& other); // Matrix multiplication assignment with transposed matrix (from the left)
	Matrix& operator/=(Scalar scalar); // Scalar division assignment
	friend double determinant<>(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m);
	friend Matrix<ScalarParam,numColumnsParam,numRowsParam> invert<>(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m);
	ComponentArray<ScalarParam,numColumnsParam> transposeMultiply(const ComponentArray<ScalarParam,numRowsParam>& ca) const // Vector multiplication (from the left)
		{
		return MO::transposeMultiply(c,ca);
		}
	friend ComponentArray<ScalarParam,numRowsParam> operator*<>(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m,const ComponentArray<ScalarParam,numColumnsParam>& ca);
	friend ComponentArray<ScalarParam,numColumnsParam> operator/<>(const ComponentArray<ScalarParam,numRowsParam>& ca,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m);
	};

/*************************************
Operations on objects of class Matrix:
*************************************/

template <class ScalarParam,int numRowsParam,int numColumnsParam>
bool operator==(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m1,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m2); // Equality operator
template <class ScalarParam,int numRowsParam,int numColumnsParam>
bool operator!=(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m1,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m2); // Inequality operator
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator+(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m1,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m2); // Addition
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator-(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m1,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m2); // Subtraction
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator*(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m,ScalarParam scalar); // Scalar multiplication (from the right)
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator*(ScalarParam scalar,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m); // Ditto (from the left)
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator/(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m,ScalarParam scalar); // Scalar division
template <class ScalarParam,int numRowsParam,int middleParam,int numColumnsParam>
Matrix<ScalarParam,numRowsParam,numColumnsParam> operator*(const Matrix<ScalarParam,numRowsParam,middleParam>& m1,const Matrix<ScalarParam,middleParam,numColumnsParam>& m2); // Matrix multiplication
template <class ScalarParam,int numRowsParam,int numColumnsParam>
Matrix<ScalarParam,numColumnsParam,numRowsParam> transpose(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m); // Returns the transpose of a matrix

template <class ScalarParam,int numRowsParam,int numColumnsParam>
inline double determinant(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m) // Returns the determinant of a square matrix
	{
	return MatrixOperations<ScalarParam,numRowsParam,numColumnsParam>::determinant(m.c);
	}

template <class ScalarParam,int numRowsParam,int numColumnsParam>
inline Matrix<ScalarParam,numColumnsParam,numRowsParam> invert(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m) // Returns the inverse of a square matrix
	{
	Matrix<ScalarParam,numColumnsParam,numRowsParam> result;
	MatrixOperations<ScalarParam,numRowsParam,numColumnsParam>::invert(m.c,result.c);
	return result;
	}

template <class ScalarParam,int numRowsParam,int numColumnsParam>
inline ComponentArray<ScalarParam,numRowsParam> operator*(const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m,const ComponentArray<ScalarParam,numColumnsParam>& ca) // Multiplies matrix and vector
	{
	return MatrixOperations<ScalarParam,numRowsParam,numColumnsParam>::multiply(m.c,ca);
	}

template <class ScalarParam,int numRowsParam,int numColumnsParam>
inline ComponentArray<ScalarParam,numColumnsParam> operator/(const ComponentArray<ScalarParam,numRowsParam>& ca,const Matrix<ScalarParam,numRowsParam,numColumnsParam>& m) // Multiplies inverted matrix and vector
	{
	return MatrixOperations<ScalarParam,numRowsParam,numColumnsParam>::divide(ca,m.c);
	}

}

#if defined(GEOMETRY_NONSTANDARD_TEMPLATES) && !defined(GEOMETRY_MATRIX_IMPLEMENTATION)
#include <Geometry/Matrix.icpp>
#endif

#endif
