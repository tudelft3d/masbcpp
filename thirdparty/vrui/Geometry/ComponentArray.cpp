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

#include <Geometry/ComponentArray.icpp>

namespace Geometry {

/************************************************************************
Force instantiation of all standard ComponentArray classes and functions:
************************************************************************/

template class ComponentArray<int,2>;
template ComponentArray<int,2>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<int,2>::ComponentArray(const ComponentArray<double,2>&);

template class ComponentArray<int,3>;
template ComponentArray<int,3>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<int,3>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<int,3>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<int,3>::ComponentArray(const ComponentArray<float,3>&);
template ComponentArray<int,3>::ComponentArray(const ComponentArray<double,3>&);

template class ComponentArray<int,4>;
template ComponentArray<int,4>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<int,3>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<float,3>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<double,3>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<float,4>&);
template ComponentArray<int,4>::ComponentArray(const ComponentArray<double,4>&);

template class ComponentArray<float,2>;
template ComponentArray<float,2>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<float,2>::ComponentArray(const ComponentArray<double,2>&);

template class ComponentArray<float,3>;
template ComponentArray<float,3>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<float,3>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<float,3>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<float,3>::ComponentArray(const ComponentArray<int,3>&);
template ComponentArray<float,3>::ComponentArray(const ComponentArray<double,3>&);

template class ComponentArray<float,4>;
template ComponentArray<float,4>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<int,3>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<float,3>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<double,3>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<int,4>&);
template ComponentArray<float,4>::ComponentArray(const ComponentArray<double,4>&);

template class ComponentArray<double,2>;
template ComponentArray<double,2>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<double,2>::ComponentArray(const ComponentArray<float,2>&);

template class ComponentArray<double,3>;
template ComponentArray<double,3>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<double,3>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<double,3>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<double,3>::ComponentArray(const ComponentArray<int,3>&);
template ComponentArray<double,3>::ComponentArray(const ComponentArray<float,3>&);

template class ComponentArray<double,4>;
template ComponentArray<double,4>::ComponentArray(const ComponentArray<int,2>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<float,2>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<double,2>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<int,3>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<float,3>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<double,3>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<int,4>&);
template ComponentArray<double,4>::ComponentArray(const ComponentArray<float,4>&);

}
