// Copyright (c) 2015
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// This file is part of masbcpp.
//
// masbcpp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// masbcpp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with masbcpp.  If not, see <http://www.gnu.org/licenses/>.

#ifndef TYPES_
#define TYPES_

#include <vrui/Geometry/Vector.h>
#include <vrui/Geometry/Point.h>

typedef float Scalar; // Scalar type for 3D points
typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
typedef std::vector<Point> PointList; // Type for 3D points
typedef Geometry::Vector<Scalar,3> Vector; // Type for 3D vectors
typedef std::vector<Vector> VectorList; // Type for 3D vectors

typedef std::vector<int> intList; // Type for 3D vectors

struct ma_result {
	Point c;
	int qidx;
};

struct ma_data {
	int m;

	PointList *coords;
	VectorList *normals;
	PointList *ma_coords_in;
	PointList *ma_coords_out;
	int *ma_qidx_in;
	int *ma_qidx_out;
	VectorList *ma_bisec_in;
	VectorList *ma_bisec_out;

	float *lfs;
	bool *mask;
};

#endif