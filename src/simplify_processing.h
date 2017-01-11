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

#ifndef SIMPLIFY_PROCESSING_
#define SIMPLIFY_PROCESSING_

#include "types.h"
#include "compute_normals_processing.h"
#include "compute_ma_processing.h"

struct simplify_parameters {
   double epsilon;
   double cellsize;
   double bisec_threshold;
   int bisec_k;
   double elevation_threshold;
   double minimum_density;
   double maximum_density;
   bool true_z_dim;
   bool only_inner;
   bool squared;
   bool compute_lfs;
};


// This version of simplify takes in an already calculated ma, etc.
void simplify_lfs(simplify_parameters &input_parameters, ma_data& madata);

// This version of simplify takes in only the coords of the original point cloud.
void simplify(normals_parameters &normals_params, 
              ma_parameters &ma_params, 
              simplify_parameters &simplify_params,
              PointCloud::Ptr coords,
              bool *mask); // mask *must* be allocated ahead of time to be an array of size "coords.size()".

#endif
