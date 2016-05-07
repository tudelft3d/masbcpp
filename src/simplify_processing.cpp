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
#include <random>

// OpenMP
#ifdef WITH_OPENMP
#include <omp.h>
#endif

// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>

#ifdef VERBOSEPRINT
#include <vrui/Misc/Timer.h>
#include <iostream>
#endif

// kdtree2
#include <kdtree2/kdtree2.hpp>

// typedefs
#include "simplify_processing.h"





//==============================
//   SIMPLIFY
//==============================



void compute_lfs(ma_data &madata, double bisec_threshold, bool only_inner = true)
{
#ifdef VERBOSEPRINT
   Misc::Timer t0;
#endif

   int N = 2 * madata.m;
   if (only_inner) {
      N = madata.m;
      (*madata.ma_coords).resize(N); // HACK this will destroy permanently the exterior ma_coords!
   }
   // compute bisector and filter .. rebuild kdtree .. compute lfs .. compute grid .. thin each cell
   bool ma_coords_mask[N];
   VectorList ma_bisec(N);
   //madata.ma_bisec = &ma_bisec;
   for (int i = 0; i < N; i++) {
      if (madata.ma_qidx[i] != -1) {
         Vector f1_in = (*madata.coords)[i%madata.m] - (*madata.ma_coords)[i];
         Vector f2_in = (*madata.coords)[madata.ma_qidx[i]] - (*madata.ma_coords)[i];

         ma_bisec[i] = (f1_in + f2_in).normalize();
      }
   }
#ifdef VERBOSEPRINT
   t0.elapse();
   std::cout << "Computed bisectors in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif

   int k = 2, count = 0;
   {
      kdtree2::KDTree kd_tree(*madata.ma_coords, true);
      kd_tree.sort_results = true;
#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Constructed kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif

      kdtree2::KDTreeResultVector result;
#pragma omp parallel for private(result)
      for (int i = 0; i < N; i++) {
         ma_coords_mask[i] = false;
         if (madata.ma_qidx[i] != -1) {
            kd_tree.n_nearest((*madata.ma_coords)[i], k, result);

            float bisec_angle = acos(ma_bisec[result[1].idx] * ma_bisec[i]);
            if (bisec_angle < bisec_threshold)
               ma_coords_mask[i] = true;
         }
      }
      for (int i = 0; i < N; i++)
         if (ma_coords_mask[i])
            count++;

#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Cleaned MA points in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif
   }
   // mask and copy pointlist ma_coords
   PointList ma_coords_masked(count);
   int j = 0;
   for (int i = 0; i < N; i++) {
      if (ma_coords_mask[i])
         ma_coords_masked[j++] = (*madata.ma_coords)[i];
   }
#ifdef VERBOSEPRINT
   t0.elapse();
   std::cout << "Copied cleaned MA points in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif

   k = 1;
   {
      // rebuild kd-tree
      kdtree2::KDTree kd_tree(ma_coords_masked, true);
      // kd_tree = new kdtree2::KDTree;
      kd_tree.sort_results = true;
#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Constructed cleaned kd-tree in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif

      kdtree2::KDTreeResultVector result;
#pragma omp parallel for private(result)
      for (int i = 0; i < madata.m; i++) {
         kd_tree.n_nearest((*madata.coords)[i], k, result);
         madata.lfs[i] = sqrt(result[0].dis);
      }
#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Computed LFS in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif
   }


}

inline int flatindex(int ind[], int size[], int dimension) {
   if (dimension == 2)
      return ind[0] + size[0] * ind[1];
   return ind[0] + size[0] * (ind[1] + ind[2] * size[1]);
}

void simplify(ma_data &madata, double cellsize, double epsilon, int dimension = 3, double elevation_threshold = 0.0) {
#ifdef VERBOSEPRINT
   Misc::Timer t0;
#endif

   Box::Size size = madata.bbox.getSize();
   Point origin = Point(madata.bbox.min);

   int* resolution = new int[dimension];

   #ifdef VERBOSEPRINT
   std::cout << "Grid resolution: ";
   #endif
   for (int i = 0; i < dimension; i++) {
      resolution[i] = int(size[i] / cellsize) + 1;
      #ifdef VERBOSEPRINT
      std::cout << resolution[i] << " ";
      #endif
   }
   #ifdef VERBOSEPRINT
   std::cout << std::endl;
   #endif

   int ncells = 1;
   for (int i = 0; i < dimension; i++)
      ncells *= resolution[i];

   intList** grid = new intList*[ncells];
   for (int i = 0; i < ncells; i++) {
      grid[i] = NULL;
   }

   int* idx = new int[dimension];
   int index;
   for (int i = 0; i < madata.m; i++) {
      for (int j = 0; j < dimension; j++) {
         idx[j] = int(((*madata.coords)[i][j] - origin[j]) / cellsize);
      }
      index = flatindex(idx, resolution, dimension);

      if (grid[index] == NULL) {
         intList *ilist = new intList;
         // std::unique_ptr<intList> ilist(new intList);
         grid[index] = ilist;
      }
      (*grid[index]).push_back(i);
   }

   delete[] resolution; resolution = NULL;
   delete[] idx; idx = NULL;

   double mean_lfs, target_n, A = cellsize*cellsize;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<float> randu(0, 1);

   // parallelize?
   for (int i = 0; i < ncells; i++)
      if (grid[i] != NULL) {
         size_t n = grid[i]->size();
         float sum = 0, max_z, min_z;
         max_z = min_z = (*madata.coords)[(*grid[i])[0]][2];

         for (auto j : *grid[i]) {
            sum += madata.lfs[j];
            float z = (*madata.coords)[j][2];
            if (z > max_z) max_z = z;
            if (z < min_z) min_z = z;
         }

         mean_lfs = sum / n;

         if (elevation_threshold != 0 && (max_z - min_z) > elevation_threshold)
            mean_lfs /= 5;


         target_n = A / pow(epsilon*mean_lfs, 2);
         for (auto j : *grid[i])
            madata.mask[j] = randu(gen) <= target_n / n;
      }
#ifdef VERBOSEPRINT
   t0.elapse();
   std::cout << "Performed grid simplification in " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif


   // clear some memory in non-smart ptr way
   for (int i = 0; i < ncells; i++) {
      delete grid[i];
   }
   delete[] grid;

}

void simplify_lfs(simplify_parameters &input_parameters, ma_data& madata)
{

   // compute lfs, simplify
   compute_lfs(madata, input_parameters.bisec_threshold, input_parameters.only_inner);
   simplify(madata, input_parameters.cellsize, input_parameters.epsilon, input_parameters.dimension, input_parameters.elevation_threshold);
}

