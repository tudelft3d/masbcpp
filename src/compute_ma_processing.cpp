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

#include "compute_ma_processing.h"

#include <limits>

#ifdef VERBOSEPRINT
#include <chrono>
#include <iostream>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#ifdef VERBOSEPRINT
typedef std::chrono::high_resolution_clock Clock;
#endif

//==============================
//   COMPUTE MA
//==============================

const Scalar delta_convergance = 1E-5f;
const unsigned int iteration_limit = 30;
const Vector3 nanPoint(std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::quiet_NaN());

inline Scalar compute_radius(const Vector3 &p, const Vector3 &n, const Vector3 &q) {
   // Compute radius of the ball that touches points p and q and whose center falls on the normal n from p
   Scalar d = (p - q).norm();
   Scalar cos_theta = n.dot(p - q) / d;
   return Scalar(d / (2 * cos_theta));
}

inline Scalar cos_angle(const Vector3 p, const Vector3 q) {
   // Calculate the cosine of angle between vector p and q, see http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
   Scalar result = p.dot(q) / (p.norm() * q.norm());
   if (result > 1) return 1;
   else if (result < -1) return -1;
   return result;
}

ma_result sb_point(const ma_parameters &input_parameters, const Vector3 &p, const Vector3 &n, pcl::search::KdTree<Point>::Ptr kd_tree) {
   // Calculate a medial ball for a given oriented point using the shrinking ball algorithm,
   // see https://3d.bk.tudelft.nl/rypeters/pdfs/16candg.pdf section 3.2 for details
   unsigned int j = 0;
   Scalar r = input_parameters.initial_radius, d;
   Vector3 q, c_next;
   int qidx = -1, qidx_next;
   Point c; c.getVector3fMap() = p - n * r;

   // Results from our search
   std::vector<int> k_indices(2);
   std::vector<Scalar> k_distances(2);

   while (true) {

      // find closest point to c
      kd_tree->nearestKSearch(c, 1, k_indices, k_distances);

      qidx_next = k_indices[0];
      q = kd_tree->getInputCloud()->at(qidx_next).getVector3fMap();
      d = k_distances[0];
      
      // This should handle all (special) cases where we want to break the loop
      // - normal case when ball no longer shrinks
      // - the case where q==p
      // - any duplicate point cases
      if ((d >= (r*r)-delta_convergance) || (p==q))
            break;

      // compute new radius
      r = compute_radius(p, n, q);

      // compute next ball center
      c_next = p - n * r;

      // denoising
      if (input_parameters.denoise_preserve || input_parameters.denoise_planar) {
         Scalar a = cos_angle(p - c_next, q - c_next);
         Scalar separation_angle = std::acos(a);

         if (input_parameters.denoise_preserve && (separation_angle < input_parameters.denoise_preserve && j>0 && r > (q - p).norm())) {
            break;
         }
         if (input_parameters.denoise_planar && (separation_angle < input_parameters.denoise_planar && j == 0)) {
            c.getVector3fMap() = input_parameters.nan_for_initr ? nanPoint : p - n * input_parameters.initial_radius;
            break;
         }
      }

      // stop iteration if this looks like an infinite loop:
      if (j > iteration_limit)
         break;

      c.getVector3fMap() = c_next;
      qidx = qidx_next;
      j++;
   }

   return{ c, qidx };
}

void sb_points(ma_parameters &input_parameters, ma_data &madata, bool inner = 1) {
   // outer mat should be written to second half of ma_coords/ma_qidx
   size_t offset = 0;
   if (inner == false)
      offset = madata.coords->size();

#pragma omp parallel for
   for (int i = 0; i < madata.coords->size(); i++) {
      Vector3 p = (*madata.coords)[i].getVector3fMap();
      Vector3 n;
      if (inner)
         n = (*madata.normals)[i].getNormalVector3fMap();
      else
         n = -(*madata.normals)[i].getNormalVector3fMap();

      ma_result r = sb_point(input_parameters, p, n, madata.kd_tree);

      (*madata.ma_coords)[i + offset] = r.c;
      madata.ma_qidx[i + offset] = r.qidx;
   }
}

void compute_masb_points(ma_parameters &input_parameters, ma_data &madata) {
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   if (!madata.kd_tree) {
      madata.kd_tree.reset(new pcl::search::KdTree<Point>());
      madata.kd_tree->setInputCloud(madata.coords);
#ifdef VERBOSEPRINT
      auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Constructed kd-tree in " << elapsed_time.count() << " ms" << std::endl;
      start_time = Clock::now();
#endif
   }

   // Inside processing
   sb_points(input_parameters, madata, 1);
#ifdef VERBOSEPRINT
   auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Done shrinking interior balls, took " << elapsed_time.count() << " ms" << std::endl;
   start_time = Clock::now();
#endif

   // Outside processing
   sb_points(input_parameters, madata, 0);
#ifdef VERBOSEPRINT
   elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Done shrinking exterior balls, took " << elapsed_time.count() << " ms" << std::endl;
#endif
}

