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

#include <pcl/common/common.h>

#ifdef VERBOSEPRINT
#include <chrono>
#include <iostream>
#endif

// OpenMP
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#ifdef VERBOSEPRINT
typedef std::chrono::high_resolution_clock Clock;
#endif

/*
// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>

#ifdef VERBOSEPRINT
#include <vrui/Misc/Timer.h>
#endif

// kdtree2
#include <kdtree2/kdtree2.hpp>
*/

// typedefs
#include "simplify_processing.h"





//==============================
//   SIMPLIFY
//==============================



bool compute_lfs(ma_data &madata, double bisec_threshold, int bisec_k, bool only_inner = true)
{
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   size_t N = 2 * madata.coords->size();
   if (only_inner) {
      N = madata.coords->size();
      (*madata.ma_coords).resize(N); // HACK this will destroy permanently the exterior ma_coords!
   }
   // compute bisector and filter .. rebuild kdtree .. compute lfs .. compute grid .. thin each cell

   Vector3List ma_bisec(N);
   //madata.ma_bisec = &ma_bisec;
   for (int i = 0; i < N; i++) {
      if (madata.ma_qidx[i] != -1) {
         Vector3 f1 = (*madata.coords)[i%madata.coords->size()].getVector3fMap() - (*madata.ma_coords)[i].getVector3fMap();
         Vector3 f2 = (*madata.coords)[madata.ma_qidx[i]].getVector3fMap() - (*madata.ma_coords)[i].getVector3fMap();

         ma_bisec[i] = (f1 + f2).normalized();
      }
   }
#ifdef VERBOSEPRINT
   auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Computed bisectors in " << elapsed_time.count() << " ms" << std::endl;
   start_time = Clock::now();
#endif

   int count = 0;
   std::vector<bool> bisec_mask(N);
   {
      pcl::search::KdTree<Point>::Ptr kd_tree(new pcl::search::KdTree<Point>());
      kd_tree->setInputCloud(madata.ma_coords);
#ifdef VERBOSEPRINT
      auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Constructed kd-tree in " << elapsed_time.count() << " ms" << std::endl;
      start_time = Clock::now();
#endif

      // Results from our search
      std::vector<int> k_indices(bisec_k);
      std::vector<Scalar> k_distances(bisec_k);

#pragma omp parallel for private(k_indices, k_distances)
      for (int i = 0; i < N; i++) {
         bisec_mask[i] = false;
         if (madata.ma_qidx[i] != -1) {
            kd_tree->nearestKSearch((*madata.ma_coords)[i], bisec_k, k_indices, k_distances); // find closest point to c

            float bisec_angle, max_bisec_angle = 0;
            for (int j=1; j<bisec_k; j++){
                  bisec_angle = std::acos(ma_bisec[k_indices[j]].dot(ma_bisec[i]));
                  if (bisec_angle > max_bisec_angle)
                        max_bisec_angle = bisec_angle;
            }
            if (max_bisec_angle < bisec_threshold)
               bisec_mask[i] = true;
         }
      }
      for (int i = 0; i < N; i++)
         if (bisec_mask[i])
            count++;

#ifdef VERBOSEPRINT
      elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Cleaned MA points in " << elapsed_time.count() << " ms" << std::endl;
      start_time = Clock::now();
#endif
   }

   // We can't produce LFS values if there are no MAT points
   if (count == 0)
      return false;

   // mask and copy pointlist ma_coords
   PointCloud::Ptr ma_coords_masked(new PointCloud);
   ma_coords_masked->reserve(count);

   for (int i = 0; i < N; i++) {
      if (bisec_mask[i])
         ma_coords_masked->push_back((*madata.ma_coords)[i]);
   }
#ifdef VERBOSEPRINT
   elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Copied cleaned MA points in " << elapsed_time.count() << " ms" << std::endl;
   start_time = Clock::now();
#endif

   {
      // rebuild kd-tree
      pcl::search::KdTree<Point>::Ptr kd_tree(new pcl::search::KdTree<Point>());
      kd_tree->setInputCloud(ma_coords_masked);
      // kd_tree = new kdtree2::KDTree;
#ifdef VERBOSEPRINT
      elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Constructed cleaned kd-tree in " << elapsed_time.count() << " ms" << std::endl;
      start_time = Clock::now();
#endif

      // Results from our search
      std::vector<int> k_indices(1);
      std::vector<Scalar> k_distances(1);

#pragma omp parallel for private(k_indices, k_distances)
      for (int i = 0; i < madata.coords->size(); i++) {
         kd_tree->nearestKSearch((*madata.coords)[i], 1, k_indices, k_distances); // find closest point to c

         madata.lfs[i] = std::sqrt(k_distances[0]);
      }
#ifdef VERBOSEPRINT
      elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Computed LFS in " << elapsed_time.count() << " ms" << std::endl;
      start_time = Clock::now();
#endif
   }

   return true;
}

inline int flatindex(int ind[], int size[], bool true_z_dim) {
   if (!true_z_dim)
      return ind[0] + size[0] * ind[1];
   return ind[0] + size[0] * (ind[1] + ind[2] * size[1]);
}

void simplify(ma_data &madata, 
             double cellsize, 
             double epsilon, 
             bool true_z_dim = true, 
             double elevation_threshold = 0.0, 
             double minimum_density = 0,
             double maximum_density = 0,
             bool squared = false) 
{
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   Point minPt;
   Point maxPt;
   pcl::getMinMax3D(*(madata.coords), minPt, maxPt);
   float size[3];
   size[0] = maxPt.x - minPt.x;
   size[1] = maxPt.y - minPt.y;
   if (true_z_dim)
      size[2] = maxPt.z - minPt.z;
   Point origin = minPt;

   //Box::Size size = madata.bbox.getSize();
   //Point origin = Point(madata.bbox.min);

   int* resolution = new int[3];

   #ifdef VERBOSEPRINT
   std::cout << "Epsilon: " << epsilon << std::endl;
   std::cout << "Maximum density: " << maximum_density << std::endl;
   std::cout << "Minimum density: " << minimum_density << std::endl;
   std::cout << "True z: " << true_z_dim << std::endl;
   std::cout << "Squared: " << squared << std::endl;
   std::cout << "Cellsize: " << cellsize << std::endl;
   std::cout << "Data dimensions: " << size[0] << " x " << size[1];
   if (true_z_dim) std::cout << " x " << size[2];
   std::cout << std::endl;
   std::cout << "Grid dimension: ";
   #endif

   // x, y, z - resolution
   resolution[0] = int(size[0] / cellsize) + 1;
   resolution[1] = int(size[1] / cellsize) + 1;
   if (true_z_dim)
      resolution[2] = int(size[2] / cellsize) + 1;

   #ifdef VERBOSEPRINT
   std::cout << resolution[0] << " ";
   std::cout << resolution[1] << " ";
   if (true_z_dim)
      std::cout << resolution[2] << " ";
   std::cout << std::endl;
   #endif

   int ncells = 1;
   ncells *= resolution[0];
   ncells *= resolution[1];
   if (true_z_dim)
      ncells *= resolution[2];

   intList** grid = new intList*[ncells];
   for (size_t i = 0; i < ncells; i++) {
      grid[i] = NULL;
   }

   int* idx = new int[3];
   int index;
   for (int i = 0; i < madata.coords->size(); i++) {
      idx[0] = int(((*madata.coords)[i].x - origin.x) / cellsize);
      idx[1] = int(((*madata.coords)[i].y - origin.y) / cellsize);
      if (true_z_dim)
         idx[2] = int(((*madata.coords)[i].z - origin.z) / cellsize);

      index = flatindex(idx, resolution, true_z_dim);

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

#ifdef DETERMINISTIC_RNG
   std::mt19937 gen;
#else
   std::random_device rd;
   std::mt19937 gen(rd());
#endif
   std::uniform_real_distribution<float> randu(0, 1);

   double target_n_max = maximum_density * A;
   double target_n_min = minimum_density * A;
   // parallelize?
   for (int i = 0; i < ncells; i++)
      if (grid[i] != NULL) {
         size_t n = grid[i]->size();
         float sum = 0, max_z, min_z;
         max_z = min_z = (*madata.coords)[(*grid[i])[0]].z;

         for (auto j : *grid[i]) {
            sum += madata.lfs[j];
            float z = (*madata.coords)[j].z;
            if (z > max_z) max_z = z;
            if (z < min_z) min_z = z;
         }

         mean_lfs = sum / n;

         if (squared) mean_lfs = pow(mean_lfs, 2);
         if (elevation_threshold != 0 && (max_z - min_z) > elevation_threshold)
            mean_lfs /= 10;
            // mean_lfs = 0.01;

         target_n = A / pow(epsilon*mean_lfs, 2);
         if(target_n_max != 0 && target_n > target_n_max) target_n = target_n_max;
         else if(target_n_min != 0 && target_n < target_n_min) target_n = target_n_min;
         for (auto j : *grid[i])
            madata.mask[j] = randu(gen) <= target_n / n;
      }
#ifdef VERBOSEPRINT
   auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Performed grid simplification in " << elapsed_time.count() << " ms" << std::endl;
   start_time = Clock::now();
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
   if (input_parameters.compute_lfs)
   {
      // If we can't compute LFS values, leave the mask as all false
      if (!compute_lfs(madata, input_parameters.bisec_threshold, input_parameters.bisec_k, input_parameters.only_inner))
         return;
   }
   simplify(madata, input_parameters.cellsize, 
                    input_parameters.epsilon, 
                    input_parameters.true_z_dim, 
                    input_parameters.elevation_threshold, 
                    input_parameters.minimum_density,
                    input_parameters.maximum_density,
                    input_parameters.squared);
}

void simplify(normals_parameters &normals_params, 
              ma_parameters &ma_params,
              simplify_parameters &simplify_params,
              PointCloud::Ptr coords, bool *mask) // mask *must* be allocated ahead of time to be an array of size "coords.size()".
{
   ///////////////////////////
   // Step 0: prepare data struct:
   ma_data madata = {};
   madata.coords = coords; // add to the reference count
   

   ///////////////////////////
   // Step 1: compute normals:
   NormalCloud::Ptr normals(new NormalCloud);
   normals->resize(madata.coords->size());
   madata.normals = normals; // add to the reference count
   compute_normals(normals_params, madata);


   ///////////////////////////
   // Step 2: compute ma
   PointCloud::Ptr ma_coords(new PointCloud);
   ma_coords->resize(2*madata.coords->size());
   madata.ma_coords = ma_coords; // add to the reference count
   madata.ma_qidx.resize(2 * madata.coords->size());
   compute_masb_points(ma_params, madata);

   //delete madata.kdtree_coords; madata.kdtree_coords = NULL;


   ///////////////////////////
   // Step 3: Simplify
   madata.mask.resize(madata.coords->size());
   madata.lfs.resize(madata.coords->size());
   simplify_lfs(simplify_params, madata);

   ///////////////////////////
   // Pass back the results in a safe way.
   std::copy(madata.mask.begin(), madata.mask.end(), mask);
   //memcpy(mask, &madata.mask[0], madata.mask.size() * sizeof(bool));
}

