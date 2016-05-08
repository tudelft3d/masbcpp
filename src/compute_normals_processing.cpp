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
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>
#include <vrui/Geometry/PCACalculator.h>

#ifdef VERBOSEPRINT
#include <vrui/Misc/Timer.h>
#include <iostream>
#endif

// typedefs
#include "compute_normals_processing.h"





//==============================
//   COMPUTE NORMALS
//==============================


// Vector estimate_normal(Point &p, kdtree* kd_tree, int k)
// {
//    kdtree_result result(k);
//    kd_tree->search(result, &p[0], 1, k+1);

//    Geometry::PCACalculator<3> PCACalc;
//    for (int i = 0; i < k + 1; i++)
//       PCACalc.accumulatePoint(Point(&(kd_tree->points[ result.idx[i]*3 ])));
//       // PCACalc.accumulatePoint(kd_tree->the_data[result[i].idx]);

//    double eigen_values[3];
//    PCACalc.calcCovariance();
//    PCACalc.calcEigenvalues(eigen_values);
//    return PCACalc.calcEigenvector(eigen_values[2]);
// }

void estimate_normals(ma_data &madata, int k)
{
      omp_set_nested(0); // disable kdtree parallelism 
      #pragma omp parallel for
      for (int i = 0; i < madata.m; i++) {
            // Vector n = estimate_normal(Point(&(madata.coords[i*3])), madata.kdtree_coords, k);
            float *p = &(madata.coords[i*3]);
  
            kdtree_result result(k);
            madata.kdtree_coords->search(result, p, 1, k+1);

            Geometry::PCACalculator<3> PCACalc;
            for (int i = 0; i < k + 1; i++)
                  PCACalc.accumulatePoint(Point(&(madata.coords[ result.idx[i]*3 ])));
                  // PCACalc.accumulatePoint(kd_tree->the_data[result[i].idx]);

            double eigen_values[3];
            PCACalc.calcCovariance();
            PCACalc.calcEigenvalues(eigen_values);
            Vector n = PCACalc.calcEigenvector(eigen_values[2]);
            madata.normals[i*3+0] = n[0];
            madata.normals[i*3+1] = n[1];
            madata.normals[i*3+2] = n[2];
      }
}

void compute_normals(normals_parameters &input_parameters, ma_data &madata)
{
      #ifdef VERBOSEPRINT
      Misc::Timer t0;
      #endif
      
      if (madata.kdtree_coords == NULL) {
            // madata.kdtree_coords = new kdtree2::KDTree((*madata.coords), input_parameters.kd_tree_reorder);
            madata.kdtree_coords = new kdtree(madata.coords, 3, madata.m, 16);
            #ifdef VERBOSEPRINT
            t0.elapse();
            std::cout << "Constructed kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
            #endif
      }
      // madata.kdtree_coords->sort_results = false;

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

