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

#ifndef MASBCPP_MADATA_
#define MASBCPP_MADATA_

// kdtree2
#include "types.h"
#include <kdtree2/kdtree2.hpp>

struct ma_data {
   unsigned int m;
//    Box bbox;

   PointCloud::Ptr coords; // don't own this memory
   NormalCloud::Ptr normals; // don't own this memory
   PointCloud::Ptr ma_coords; // don't own this memory
   int ma_qidx;

   float lfs;
   bool mask;
   
   kdtree2::KDTree* kdtree_coords;
};

struct io_parameters {
    bool coords;
    bool normals;
    bool ma_coords;
    bool ma_qidx;
};

#endif