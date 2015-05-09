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

#include <Geometry/Matrix.icpp>

namespace Geometry {

/****************************************************************
Force instantiation of all standard Matrix classes and functions:
****************************************************************/

template class MatrixOperations<float,2,2>;
template class MatrixOperations<float,2,3>;
template class MatrixOperations<float,3,3>;
template class MatrixOperations<float,3,4>;
template class MatrixOperations<float,4,4>;
template class MatrixOperations<double,2,2>;
template class MatrixOperations<double,2,3>;
template class MatrixOperations<double,3,3>;
template class MatrixOperations<double,3,4>;
template class MatrixOperations<double,4,4>;

template class Matrix<float,2,2>;
template Matrix<float,2,2> Matrix<float,2,2>::fromRowMajor(const float*);
template Matrix<float,2,2> Matrix<float,2,2>::fromRowMajor(const double*);
template Matrix<float,2,2> Matrix<float,2,2>::fromColumnMajor(const float*);
template Matrix<float,2,2> Matrix<float,2,2>::fromColumnMajor(const double*);
template Matrix<float,2,2>::Matrix(const Matrix<double,2,2>&);
template Matrix<float,2,2> operator+(const Matrix<float,2,2>&,const Matrix<float,2,2>&);
template Matrix<float,2,2> operator-(const Matrix<float,2,2>&,const Matrix<float,2,2>&);
template Matrix<float,2,2> operator*(const Matrix<float,2,2>&,float);
template Matrix<float,2,2> operator*(float,const Matrix<float,2,2>&);
template Matrix<float,2,2> operator/(const Matrix<float,2,2>&,float);
template Matrix<float,2,2> operator*(const Matrix<float,2,2>&,const Matrix<float,2,2>&);
template Matrix<float,2,2> transpose(const Matrix<float,2,2>&);

template class Matrix<float,2,3>;
template Matrix<float,2,3> Matrix<float,2,3>::fromRowMajor(const float*);
template Matrix<float,2,3> Matrix<float,2,3>::fromRowMajor(const double*);
template Matrix<float,2,3> Matrix<float,2,3>::fromColumnMajor(const float*);
template Matrix<float,2,3> Matrix<float,2,3>::fromColumnMajor(const double*);
template Matrix<float,2,3>::Matrix(const Matrix<double,2,3>&);
template Matrix<float,2,3> operator+(const Matrix<float,2,3>&,const Matrix<float,2,3>&);
template Matrix<float,2,3> operator-(const Matrix<float,2,3>&,const Matrix<float,2,3>&);
template Matrix<float,2,3> operator*(const Matrix<float,2,3>&,float);
template Matrix<float,2,3> operator*(float,const Matrix<float,2,3>&);
template Matrix<float,2,3> operator/(const Matrix<float,2,3>&,float);
template Matrix<float,2,3> operator*(const Matrix<float,2,3>&,const Matrix<float,3,3>&);
template Matrix<float,2,3> operator*(const Matrix<float,2,2>&,const Matrix<float,2,3>&);

template class Matrix<float,3,3>;
template Matrix<float,3,3> Matrix<float,3,3>::fromRowMajor(const float*);
template Matrix<float,3,3> Matrix<float,3,3>::fromRowMajor(const double*);
template Matrix<float,3,3> Matrix<float,3,3>::fromColumnMajor(const float*);
template Matrix<float,3,3> Matrix<float,3,3>::fromColumnMajor(const double*);
template Matrix<float,3,3>::Matrix(const Matrix<double,3,3>&);
template Matrix<float,3,3> operator+(const Matrix<float,3,3>&,const Matrix<float,3,3>&);
template Matrix<float,3,3> operator-(const Matrix<float,3,3>&,const Matrix<float,3,3>&);
template Matrix<float,3,3> operator*(const Matrix<float,3,3>&,float);
template Matrix<float,3,3> operator*(float,const Matrix<float,3,3>&);
template Matrix<float,3,3> operator/(const Matrix<float,3,3>&,float);
template Matrix<float,3,3> operator*(const Matrix<float,3,3>&,const Matrix<float,3,3>&);
template Matrix<float,3,3> transpose(const Matrix<float,3,3>&);

template class Matrix<float,3,4>;
template Matrix<float,3,4> Matrix<float,3,4>::fromRowMajor(const float*);
template Matrix<float,3,4> Matrix<float,3,4>::fromRowMajor(const double*);
template Matrix<float,3,4> Matrix<float,3,4>::fromColumnMajor(const float*);
template Matrix<float,3,4> Matrix<float,3,4>::fromColumnMajor(const double*);
template Matrix<float,3,4>::Matrix(const Matrix<double,3,4>&);
template Matrix<float,3,4> operator+(const Matrix<float,3,4>&,const Matrix<float,3,4>&);
template Matrix<float,3,4> operator-(const Matrix<float,3,4>&,const Matrix<float,3,4>&);
template Matrix<float,3,4> operator*(const Matrix<float,3,4>&,float);
template Matrix<float,3,4> operator*(float,const Matrix<float,3,4>&);
template Matrix<float,3,4> operator/(const Matrix<float,3,4>&,float);
template Matrix<float,3,4> operator*(const Matrix<float,3,4>&,const Matrix<float,4,4>&);
template Matrix<float,3,4> operator*(const Matrix<float,3,3>&,const Matrix<float,3,4>&);

template class Matrix<float,4,4>;
template Matrix<float,4,4> Matrix<float,4,4>::fromRowMajor(const float*);
template Matrix<float,4,4> Matrix<float,4,4>::fromRowMajor(const double*);
template Matrix<float,4,4> Matrix<float,4,4>::fromColumnMajor(const float*);
template Matrix<float,4,4> Matrix<float,4,4>::fromColumnMajor(const double*);
template Matrix<float,4,4>::Matrix(const Matrix<double,4,4>&);
template Matrix<float,4,4> operator+(const Matrix<float,4,4>&,const Matrix<float,4,4>&);
template Matrix<float,4,4> operator-(const Matrix<float,4,4>&,const Matrix<float,4,4>&);
template Matrix<float,4,4> operator*(const Matrix<float,4,4>&,float);
template Matrix<float,4,4> operator*(float,const Matrix<float,4,4>&);
template Matrix<float,4,4> operator/(const Matrix<float,4,4>&,float);
template Matrix<float,4,4> operator*(const Matrix<float,4,4>&,const Matrix<float,4,4>&);
template Matrix<float,4,4> transpose(const Matrix<float,4,4>&);

template class Matrix<double,2,2>;
template Matrix<double,2,2> Matrix<double,2,2>::fromRowMajor(const float*);
template Matrix<double,2,2> Matrix<double,2,2>::fromRowMajor(const double*);
template Matrix<double,2,2> Matrix<double,2,2>::fromColumnMajor(const float*);
template Matrix<double,2,2> Matrix<double,2,2>::fromColumnMajor(const double*);
template Matrix<double,2,2>::Matrix(const Matrix<float,2,2>&);
template Matrix<double,2,2> operator+(const Matrix<double,2,2>&,const Matrix<double,2,2>&);
template Matrix<double,2,2> operator-(const Matrix<double,2,2>&,const Matrix<double,2,2>&);
template Matrix<double,2,2> operator*(const Matrix<double,2,2>&,double);
template Matrix<double,2,2> operator*(double,const Matrix<double,2,2>&);
template Matrix<double,2,2> operator/(const Matrix<double,2,2>&,double);
template Matrix<double,2,2> operator*(const Matrix<double,2,2>&,const Matrix<double,2,2>&);
template Matrix<double,2,2> transpose(const Matrix<double,2,2>&);

template class Matrix<double,2,3>;
template Matrix<double,2,3> Matrix<double,2,3>::fromRowMajor(const float*);
template Matrix<double,2,3> Matrix<double,2,3>::fromRowMajor(const double*);
template Matrix<double,2,3> Matrix<double,2,3>::fromColumnMajor(const float*);
template Matrix<double,2,3> Matrix<double,2,3>::fromColumnMajor(const double*);
template Matrix<double,2,3>::Matrix(const Matrix<float,2,3>&);
template Matrix<double,2,3> operator+(const Matrix<double,2,3>&,const Matrix<double,2,3>&);
template Matrix<double,2,3> operator-(const Matrix<double,2,3>&,const Matrix<double,2,3>&);
template Matrix<double,2,3> operator*(const Matrix<double,2,3>&,double);
template Matrix<double,2,3> operator*(double,const Matrix<double,2,3>&);
template Matrix<double,2,3> operator/(const Matrix<double,2,3>&,double);
template Matrix<double,2,3> operator*(const Matrix<double,2,3>&,const Matrix<double,3,3>&);
template Matrix<double,2,3> operator*(const Matrix<double,2,2>&,const Matrix<double,2,3>&);

template class Matrix<double,3,3>;
template Matrix<double,3,3> Matrix<double,3,3>::fromRowMajor(const float*);
template Matrix<double,3,3> Matrix<double,3,3>::fromRowMajor(const double*);
template Matrix<double,3,3> Matrix<double,3,3>::fromColumnMajor(const float*);
template Matrix<double,3,3> Matrix<double,3,3>::fromColumnMajor(const double*);
template Matrix<double,3,3>::Matrix(const Matrix<float,3,3>&);
template Matrix<double,3,3> operator+(const Matrix<double,3,3>&,const Matrix<double,3,3>&);
template Matrix<double,3,3> operator-(const Matrix<double,3,3>&,const Matrix<double,3,3>&);
template Matrix<double,3,3> operator*(const Matrix<double,3,3>&,double);
template Matrix<double,3,3> operator*(double,const Matrix<double,3,3>&);
template Matrix<double,3,3> operator/(const Matrix<double,3,3>&,double);
template Matrix<double,3,3> operator*(const Matrix<double,3,3>&,const Matrix<double,3,3>&);
template Matrix<double,3,3> transpose(const Matrix<double,3,3>&);

template class Matrix<double,3,4>;
template Matrix<double,3,4> Matrix<double,3,4>::fromRowMajor(const float*);
template Matrix<double,3,4> Matrix<double,3,4>::fromRowMajor(const double*);
template Matrix<double,3,4> Matrix<double,3,4>::fromColumnMajor(const float*);
template Matrix<double,3,4> Matrix<double,3,4>::fromColumnMajor(const double*);
template Matrix<double,3,4>::Matrix(const Matrix<float,3,4>&);
template Matrix<double,3,4> operator+(const Matrix<double,3,4>&,const Matrix<double,3,4>&);
template Matrix<double,3,4> operator-(const Matrix<double,3,4>&,const Matrix<double,3,4>&);
template Matrix<double,3,4> operator*(const Matrix<double,3,4>&,double);
template Matrix<double,3,4> operator*(double,const Matrix<double,3,4>&);
template Matrix<double,3,4> operator/(const Matrix<double,3,4>&,double);
template Matrix<double,3,4> operator*(const Matrix<double,3,4>&,const Matrix<double,4,4>&);
template Matrix<double,3,4> operator*(const Matrix<double,3,3>&,const Matrix<double,3,4>&);

template class Matrix<double,4,4>;
template Matrix<double,4,4> Matrix<double,4,4>::fromRowMajor(const float*);
template Matrix<double,4,4> Matrix<double,4,4>::fromRowMajor(const double*);
template Matrix<double,4,4> Matrix<double,4,4>::fromColumnMajor(const float*);
template Matrix<double,4,4> Matrix<double,4,4>::fromColumnMajor(const double*);
template Matrix<double,4,4>::Matrix(const Matrix<float,4,4>&);
template Matrix<double,4,4> operator+(const Matrix<double,4,4>&,const Matrix<double,4,4>&);
template Matrix<double,4,4> operator-(const Matrix<double,4,4>&,const Matrix<double,4,4>&);
template Matrix<double,4,4> operator*(const Matrix<double,4,4>&,double);
template Matrix<double,4,4> operator*(double,const Matrix<double,4,4>&);
template Matrix<double,4,4> operator/(const Matrix<double,4,4>&,double);
template Matrix<double,4,4> operator*(const Matrix<double,4,4>&,const Matrix<double,4,4>&);
template Matrix<double,4,4> transpose(const Matrix<double,4,4>&);

}
