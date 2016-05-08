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

#ifndef MASBCPP_KDTREE_
#define MASBCPP_KDTREE_

#include <cfloat>
#include "_kdtree_core.h"

struct kdtree_result {
    uint32_t *idx;
    float *dists;
    unsigned int k;
    
    kdtree_result(unsigned int k): k(k) {
        idx = new uint32_t[k];
        dists = new float[k];
    };
};

class kdtree {
    public:
    kdtree(float *pa, 
           int8_t no_dims, 
           uint32_t n, 
           uint32_t bsp);
    void search(kdtree_result &result,
                float *query_points, 
                uint32_t num_points, 
                uint32_t k,  
                double distance_upper_bound=FLT_MAX, 
                double eps=0);
    float *points;
    
    private:
    int8_t no_dims;
    uint32_t n;
    Tree_float* tree;
};

#endif