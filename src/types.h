/*
Copyright (c) 2016 Ravi Peters

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MASBCPP_TYPES_
#define MASBCPP_TYPES_

#include <vector>

#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d_omp.h>

typedef float Scalar;
typedef Eigen::Matrix<Scalar, 1, 3> Vector3; // Type for 3D float points
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> ArrayX3; // Type for 3D float arrays
typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Array3X;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> ArrayXi; // Type for 1D int arrays
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> ArrayXb; // Type for 1D bool arrays

typedef pcl::PointXYZ Point;
typedef pcl::Normal Normal;
typedef pcl::PointCloud<Point> PointCloud;
typedef pcl::PointCloud<Normal> NormalCloud;

#endif
