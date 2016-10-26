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
   Scalar r, r_previous = 0;
   Vector3 q, c_next;
   int qidx = -1, qidx_next;
   Point c; c.getVector3fMap() = p - n * input_parameters.initial_radius;

   // Results from our search
   std::vector<int> k_indices(2);
   std::vector<Scalar> k_distances(2);

   while (true) {
      /*
#ifdef VERBOSEPRINT
      std::cout << "\nloop iteration: " << j << ", p = (" << p[0] << "," << p[1] << "," << p[2] << ", n = (" << n[0] << "," << n[1] << "," << n[2] << ") \n";
      std::cout << "c = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
#endif
      */

      // find closest point to c
      kd_tree->nearestKSearch(c, 2, k_indices, k_distances);

      qidx_next = k_indices[0];
      q = kd_tree->getInputCloud()->at(qidx_next).getVector3fMap();

      /*
#ifdef VERBOSEPRINT
      std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";
#endif
      */

      // handle case when q==p
      if (q == p) {
         // 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
         if (r_previous == input_parameters.initial_radius) {
            r = input_parameters.initial_radius;
            c.getVector3fMap() = input_parameters.nan_for_initr ? nanPoint : p - n * r;
            break;
         }
         // 2) otherwise just pick the second closest point
         else {
            qidx_next = k_indices[1];
            q = kd_tree->getInputCloud()->at(qidx_next).getVector3fMap();
         }
      }

      // compute radius
      if (q != p) {
         // If we only have duplicate input points to match to here, our calculations will produce NAN values and
         // we get into trouble.  If we just ignore this, we eventually jump out of this loop
         // in a healthy way.
         r = compute_radius(p, n, q);
      }

      /*
#ifdef VERBOSEPRINT
      std::cout << "r = " << r << "\n";
#endif
      */

      // if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
      if (r < 0)
         r = input_parameters.initial_radius;
      // if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
      else if (r > input_parameters.initial_radius) {
         r = input_parameters.initial_radius;
         c.getVector3fMap() = input_parameters.nan_for_initr ? nanPoint : p - n * r;
         break;
      }

      // compute next ball center
      c_next = p - n * r;

      // denoising
      if (input_parameters.denoise_preserve || input_parameters.denoise_planar) {
         Scalar a = cos_angle(p - c_next, q - c_next);
         Scalar separation_angle = std::acos(a);

         if (input_parameters.denoise_preserve && (separation_angle < input_parameters.denoise_preserve && j>0 && r > (q - p).norm())) {
            // keep previous radius:
            r = r_previous;
            // qidx = qidx_next;
            break;
         }
         if (input_parameters.denoise_planar && (separation_angle < input_parameters.denoise_planar && j == 0)) {
            r = input_parameters.initial_radius;
            c.getVector3fMap() = input_parameters.nan_for_initr ? nanPoint : p - n * r;
            // qidx = qidx_next;
            break;
         }
      }

      // stop iteration if r has converged
      if (std::abs(r_previous - r) < delta_convergance)
         break;

      // stop iteration if this looks like an infinite loop:
      if (j > iteration_limit)
         break;

      r_previous = r;
      c.getVector3fMap() = c_next;
      qidx = qidx_next;
      j++;
   }

   return{ c, qidx };
}

void sb_points(ma_parameters &input_parameters, ma_data &madata, bool inner = 1) {
   // outer mat should be written to second half of ma_coords/ma_qidx
   unsigned int offset = 0;
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

