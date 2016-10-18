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

#include <limits>

// OpenMP
#ifdef WITH_OPENMP
#include <omp.h>
#endif

// Vrui
// #include <vrui/Geometry/ComponentArray.h>
// #include <vrui/Math/Math.h>
// #include <vrui/Geometry/PCACalculator.h>

#ifdef VERBOSEPRINT
#include <vrui/Misc/Timer.h>
#include <iostream>
#endif

// typedefs
#include "compute_normals_processing.h"

#include <pcl/features/normal_3d_omp.h>

//==============================
//   COMPUTE NORMALS
//==============================


void estimate_normals(ma_data &madata, int k)
{
    // #pragma omp parallel for

    // pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);

    // for (int i = 0; i < madata.coords.rows(); i++)
    //     cloud->push_back(pcl::PointXYZ(madata.coords(i,0), madata.coords(i,1), madata.coords(i,2)));

    // Create the normal estimation class, and pass the input dataset to it
    pcl::NormalEstimationOMP<Point, Normal> ne;
    ne.setInputCloud (madata.coords);

    // Create an empty kdtree representation, and pass it to the normal estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<Point>::Ptr tree (new pcl::search::KdTree<Point> ());
    ne.setSearchMethod (tree);

    // Output datasets
    // pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);

    // Use all neighbors in a sphere of radius 3cm
    ne.setKSearch(k+1);

    // Compute the features
    ne.compute (*madata.normals);
}

void compute_normals(normals_parameters &input_parameters, ma_data &madata)
{
#ifdef VERBOSEPRINT
    Misc::Timer t0;
#endif
    
//     if (madata.kdtree_coords == NULL) {
//         madata.kdtree_coords = new kdtree2::KDTree(madata.coords, input_parameters.kd_tree_reorder);
// #ifdef VERBOSEPRINT
//         t0.elapse();
//         std::cout << "Constructed kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
// #endif
//     }
//     madata.kdtree_coords->sort_results = false;
    
    {
        estimate_normals(madata, input_parameters.k);
#ifdef VERBOSEPRINT
        t0.elapse();
        std::cout << "Done estimating normals, took " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif
    }
    
    // Free memory
    // delete madata.kdtree_coords; madata.kdtree_coords = NULL;
}

