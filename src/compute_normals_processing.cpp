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


normal_result estimate_normal(Point &p, kdtree2::KDTree* kd_tree, int k, double planefit_thres)
{
   kdtree2::KDTreeResultVector result;
   kd_tree->n_nearest(p, k + 1, result);

   Geometry::PCACalculator<3> PCACalc;
   Vector neighbor_mean = Vector(0);
   for (int i = 1; i < k + 1; i++){
      Point q = kd_tree->the_data[result[i].idx];
      PCACalc.accumulatePoint(q);
      neighbor_mean[0] = neighbor_mean[0] + q[0];
      neighbor_mean[1] = neighbor_mean[1] + q[1];
      neighbor_mean[2] = neighbor_mean[2] + q[2];
   }
   neighbor_mean = neighbor_mean/Scalar(k);

   double eigen_values[3];
   PCACalc.calcCovariance();
   PCACalc.calcEigenvalues(eigen_values);
   Vector n = PCACalc.calcEigenvector(eigen_values[2]); //it is normalised

   double d, d_mean = 0; 
   for (int i = 1; i < k + 1; i++){
      Point q = kd_tree->the_data[result[i].idx];
      d = fabs(n * (q-neighbor_mean));
      d_mean += d;
   }
   d_mean /= k;
   
   d = fabs(n * (p-neighbor_mean));
//    std::cout << (d / (d+d_mean) ) << std::endl;
   return { n, (d / (d+d_mean) > planefit_thres) };
}

void estimate_normals(ma_data &madata, int k, double planefit_thres)
{
      #pragma omp parallel for
      for (int i = 0; i < madata.coords->size(); i++){
            normal_result r = estimate_normal((*madata.coords)[i], madata.kdtree_coords, k, planefit_thres);
            (*madata.normals)[i] = r.n;
            (madata.is_outlier)[i] = r.is_outlier;
      }
}

void compute_normals(normals_parameters &input_parameters, ma_data &madata)
{
      #ifdef VERBOSEPRINT
      Misc::Timer t0;
      #endif
      
      if (madata.kdtree_coords == NULL) {
            madata.kdtree_coords = new kdtree2::KDTree((*madata.coords), input_parameters.kd_tree_reorder);
            #ifdef VERBOSEPRINT
            t0.elapse();
            std::cout << "Constructed kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
            #endif
      }
      madata.kdtree_coords->sort_results = false;

      {
            // estimate_normals(madata, input_parameters.k, input_parameters.planefit_thres);
            estimate_normals(madata, input_parameters.k, 0.8);
            #ifdef VERBOSEPRINT
            t0.elapse();
            std::cout << "Done estimating normals, took " << t0.getTime()*1000.0 << " ms" << std::endl;
            #endif
      }
      
      // Free memory
      delete madata.kdtree_coords; madata.kdtree_coords = NULL;
}

