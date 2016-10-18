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

#ifndef COMPUTE_NORMALS_PROCESSING_
#define COMPUTE_NORMALS_PROCESSING_

#include "madata.h"

struct normals_parameters {
   int k;
   bool kd_tree_reorder;
};

void compute_normals(normals_parameters &input_parameters, ma_data &madata);

// based on pcl https://github.com/PointCloudLibrary/pcl/blob/46cb8fe5589e88e36d79f9b8b8e5f4ff4fceb5de/common/include/pcl/common/impl/centroid.hpp#L488
void computeMeanAndCovarianceMatrix (ArrayX3 &cloud,
                                     Eigen::Matrix3f &covariance_matrix)
{
  // create the buffer on the stack which is much faster than using cloud[indices[i]] and centroid as a buffer
  Eigen::Matrix<Scalar, 1, 9, Eigen::RowMajor> accu = Eigen::Matrix<Scalar, 1, 9, Eigen::RowMajor>::Zero ();
  size_t point_count;

    point_count = cloud.rows ();
    // For each point in the cloud
    for (size_t i = 0; i < point_count; ++i)
    {
        accu [0] += cloud(i,0) * cloud(i,0);
        accu [1] += cloud(i,0) * cloud(i,1);
        accu [2] += cloud(i,0) * cloud(i,2);
        accu [3] += cloud(i,1) * cloud(i,1); // 4
        accu [4] += cloud(i,1) * cloud(i,2); // 5
        accu [5] += cloud(i,2) * cloud(i,2); // 8
        accu [6] += cloud(i,0);
        accu [7] += cloud(i,1);
        accu [8] += cloud(i,2);
    }

    accu /= static_cast<Scalar> (point_count);

    covariance_matrix.coeffRef (0) = accu [0] - accu [6] * accu [6];
    covariance_matrix.coeffRef (1) = accu [1] - accu [6] * accu [7];
    covariance_matrix.coeffRef (2) = accu [2] - accu [6] * accu [8];
    covariance_matrix.coeffRef (4) = accu [3] - accu [7] * accu [7];
    covariance_matrix.coeffRef (5) = accu [4] - accu [7] * accu [8];
    covariance_matrix.coeffRef (8) = accu [5] - accu [8] * accu [8];
    covariance_matrix.coeffRef (3) = covariance_matrix.coeff (1);
    covariance_matrix.coeffRef (6) = covariance_matrix.coeff (2);
    covariance_matrix.coeffRef (7) = covariance_matrix.coeff (5);
}

// based on pcl https://github.com/PointCloudLibrary/pcl/blob/46cb8fe5589e88e36d79f9b8b8e5f4ff4fceb5de/features/include/pcl/features/impl/feature.hpp
inline void
solvePlaneParameters (const Eigen::Matrix3f &covariance_matrix,
                           Scalar &nx, Scalar &ny, Scalar &nz)
{
  // Extract the smallest eigenvalue and its eigenvector
  EIGEN_ALIGN16 Eigen::Vector3f::Scalar eigen_value;
  EIGEN_ALIGN16 Eigen::Vector3f eigen_vector;
  pcl::eigen33 (covariance_matrix, eigen_value, eigen_vector);

  nx = eigen_vector [0];
  ny = eigen_vector [1];
  nz = eigen_vector [2];
}

#endif
