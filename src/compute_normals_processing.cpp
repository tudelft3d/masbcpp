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

#include "compute_normals_processing.h"

#include <limits>

#ifdef VERBOSEPRINT
#include <chrono>
#include <iostream>
#endif

#include <pcl/features/normal_3d_omp.h>

#ifdef VERBOSEPRINT
typedef std::chrono::high_resolution_clock Clock;
#endif

//==============================
//   COMPUTE NORMALS
//==============================

void estimate_normals(ma_data &madata, int k) {
   // Create the normal estimation class, and pass the input dataset to it
   pcl::NormalEstimationOMP<Point, Normal> estimation;
   estimation.setInputCloud(madata.coords);

   // Create an empty kdtree representation, and pass it to the normal estimation object.
   // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
   pcl::search::KdTree<Point>::Ptr tree(new pcl::search::KdTree<Point>());
   estimation.setSearchMethod(tree);

   // Use all neighbors in a sphere of radius 3cm
   estimation.setKSearch(k + 1);

   // Compute the features
   estimation.compute(*madata.normals);
}

void compute_normals(normals_parameters &input_parameters, ma_data &madata) {
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   estimate_normals(madata, input_parameters.k);

#ifdef VERBOSEPRINT
   auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Done estimating normals, took " << elapsed_time.count() << " ms" << std::endl;
#endif
}
