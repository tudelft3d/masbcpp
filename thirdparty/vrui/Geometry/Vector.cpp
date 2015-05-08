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

#include <Geometry/Vector.icpp>

namespace Geometry {

/****************************************************************
Force instantiation of all standard Vector classes and functions:
****************************************************************/

template const Vector<int,2> Vector<int,2>::zero;

template const Vector<int,3> Vector<int,3>::zero;

template const Vector<int,4> Vector<int,4>::zero;

template const Vector<float,2> Vector<float,2>::zero;

template const Vector<float,3> Vector<float,3>::zero;

template const Vector<float,4> Vector<float,4>::zero;

template const Vector<double,2> Vector<double,2>::zero;

template const Vector<double,3> Vector<double,3>::zero;

template const Vector<double,4> Vector<double,4>::zero;

}
