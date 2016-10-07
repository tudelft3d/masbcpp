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

// Eigen
#include <Eigen/SVD>


//==============================
//   COMPUTE NORMALS
//==============================


Vector3 estimate_normal(Vector3 p, kdtree2::KDTree* kd_tree, int k)
{
   kdtree2::KDTreeResultVector result;
   kd_tree->n_nearest(p, k+1, result);

   ArrayX3 nn_result(k+1,3);
   for (int i = 0; i < k + 1; i++)
      nn_result.row(i) = kd_tree->the_data.row(result[i].idx);
   std::cout << "obtained nn result: " << std::endl << nn_result << std::endl;

   nn_result.rowwise() -= nn_result.colwise().mean();
   std::cout << "centered nn result: " << std::endl << nn_result << std::endl;

   Eigen::JacobiSVD< ArrayX3 > svd(nn_result, Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout << "SVD:" << svd.matrixU() << std::endl;
   Vector3 n = svd.matrixU().rightCols(1);

   return n;
}

void estimate_normals(ma_data &madata, int k)
{
      // #pragma omp parallel for
      for (int i = 0; i < madata.coords.rows(); i++){
            std::cout << "Computing normal for point " << i << std::endl;
            Vector3 p = madata.coords.row(i);
            Vector3 result = estimate_normal(p, madata.kdtree_coords, k);
//          madata.normals.row(i) = result;
      }
}

void compute_normals(normals_parameters &input_parameters, ma_data &madata)
{
      #ifdef VERBOSEPRINT
      Misc::Timer t0;
      #endif
      
      if (madata.kdtree_coords == NULL) {
            madata.kdtree_coords = new kdtree2::KDTree(madata.coords, input_parameters.kd_tree_reorder);
            #ifdef VERBOSEPRINT
            t0.elapse();
            std::cout << "Constructed kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
            #endif
      }
      madata.kdtree_coords->sort_results = false;

      {
            estimate_normals(madata, input_parameters.k);
            #ifdef VERBOSEPRINT
            t0.elapse();
            std::cout << "Done estimating normals, took " << t0.getTime()*1000.0 << " ms" << std::endl;
            #endif
      }
      
      // Free memory
      delete madata.kdtree_coords; madata.kdtree_coords = NULL;
}

