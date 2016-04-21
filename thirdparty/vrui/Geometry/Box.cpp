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

#include <Geometry/Box.icpp>

namespace Geometry {

/*************************************************************
Force instantiation of all standard Box classes and functions:
*************************************************************/

template class Box<int,2>;
template Box<int,2> add(const Box<int,2>&,const Box<int,2>&);
template Box<int,2> intersect(const Box<int,2>&,const Box<int,2>&);
template class Box<int,3>;
template Box<int,3> add(const Box<int,3>&,const Box<int,3>&);
template Box<int,3> intersect(const Box<int,3>&,const Box<int,3>&);

template class Box<float,2>;
template Box<float,2> add(const Box<float,2>&,const Box<float,2>&);
template Box<float,2> intersect(const Box<float,2>&,const Box<float,2>&);
template class Box<float,3>;
template Box<float,3> add(const Box<float,3>&,const Box<float,3>&);
template Box<float,3> intersect(const Box<float,3>&,const Box<float,3>&);

template class Box<double,2>;
template Box<double,2> add(const Box<double,2>&,const Box<double,2>&);
template Box<double,2> intersect(const Box<double,2>&,const Box<double,2>&);
template class Box<double,3>;
template Box<double,3> add(const Box<double,3>&,const Box<double,3>&);
template Box<double,3> intersect(const Box<double,3>&,const Box<double,3>&);

}
