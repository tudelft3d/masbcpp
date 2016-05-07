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

#ifdef VERBOSEPRINT
#include <vrui/Misc/Timer.h>
#include <iostream>
#endif

// kdtree2
#include <kdtree2/kdtree2.hpp>

// typedefs
#include "compute_ma_processing.h"




//==============================
//   COMPUTE MA
//==============================

const Scalar delta_convergance = 1E-5f;
const unsigned int iteration_limit = 30;
const Point nanPoint(std::numeric_limits<Scalar>::quiet_NaN());


inline Scalar compute_radius(Point &p, Vector &n, Point &q)
{
   // this is basic goniometry
   double d = Geometry::mag(p - q);
   Scalar cos_theta = float((n * (p - q)) / d);
   return float(d / (2 * cos_theta));
}

inline Scalar cos_angle(Vector p, Vector q)
{
   // Calculate the cosine of angle between vector p and q, see http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
   Scalar result = float(p*q / (Geometry::mag(p) * Geometry::mag(q)));
   if (result > 1) return 1;
   else if (result < -1) return -1;
   return result;
}

ma_result sb_point(ma_parameters &input_parameters, Point &p, Vector &n, kdtree2::KDTree* kd_tree)
{
   unsigned int j = 0;
   Scalar r, r_previous = 0;
   Point q, c_next;
   int qidx = -1, qidx_next;
   Point c = p - n * input_parameters.initial_radius;

   while (1)
   {
// #ifdef VERBOSEPRINT
//       std::cout << "\nloop iteration: " << j << ", p = (" << p[0] << "," << p[1] << "," << p[2] << ", n = (" << n[0] << "," << n[1] << "," << n[2] << ") \n";

//       std::cout << "c = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
// #endif

      // find closest point to c
      kdtree2::KDTreeResultVector result;
      kd_tree->n_nearest(c, 2, result);

      qidx_next = result[0].idx;
      q = kd_tree->the_data[qidx_next];

// #ifdef VERBOSEPRINT
//       std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";
// #endif

      // handle case when q==p
      if (q == p)
      {
         // 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
         if (r_previous == input_parameters.initial_radius)
         {
            r = input_parameters.initial_radius;
            c = input_parameters.nan_for_initr ? nanPoint : p - n * r;
            break;
            // 2) otherwise just pick the second closest point
         }
         else {
            qidx_next = result[1].idx;
            q = kd_tree->the_data[qidx_next];
         }
      }

      // compute radius
      r = compute_radius(p, n, q);

// #ifdef VERBOSEPRINT
//       std::cout << "r = " << r << "\n";
// #endif

      // if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
      if (r < 0)
         r = input_parameters.initial_radius;
      // if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
      else if (r > input_parameters.initial_radius)
      {
         r = input_parameters.initial_radius;
         c = input_parameters.nan_for_initr ? nanPoint : p - n * r;
         break;
      }

      // compute next ball center
      c_next = p - n * r;

      // denoising
      if (input_parameters.denoise_preserve || input_parameters.denoise_planar)
      {
         Scalar a = cos_angle(p - c_next, q - c_next);
         Scalar separation_angle = Math::acos(a);

         if (input_parameters.denoise_preserve && (separation_angle < input_parameters.denoise_preserve && j>0 && r > Geometry::mag(q - p)))
         {
            // keep previous radius:
            r = r_previous;
            // qidx = qidx_next;
            break;
         }
         if (input_parameters.denoise_planar && (separation_angle < input_parameters.denoise_planar && j == 0))
         {
            r = input_parameters.initial_radius;
            c = input_parameters.nan_for_initr ? nanPoint : p - n * r;
            // qidx = qidx_next;
            break;
         }
      }

      // stop iteration if r has converged
      if (Math::abs(r_previous - r) < delta_convergance)
         break;

      // stop iteration if this looks like an infinite loop:
      if (j > iteration_limit)
         break;

      r_previous = r;
      c = c_next;
      qidx = qidx_next;
      j++;
   }

   return{ c, qidx };
}

void sb_points(ma_parameters &input_parameters, ma_data &madata, bool inner = 1)
{
   Point p;
   Vector n;

   // outer mat should be written to second half of ma_coords/ma_qidx
   unsigned int offset = 0;
   if (inner == false)
      offset = madata.m;
      
#pragma omp parallel for private(p, n)
   for (int i = 0; i < madata.coords->size(); i++)
   {
      p = (*madata.coords)[i];
      if (inner)
         n = (*madata.normals)[i];
      else
         n = -(*madata.normals)[i];
      ma_result r = sb_point(input_parameters, p, n, madata.kdtree_coords);
      (*madata.ma_coords)[i+offset] = r.c;
      madata.ma_qidx[i+offset] = r.qidx;
   }
   // return ma_coords;
}
//PointList &coords, VectorList &normals,PointList &ma_coords_in, int* ma_qidx_in, PointList &ma_coords_out, int* ma_qidx_out
void compute_masb_points(ma_parameters &input_parameters, ma_data &madata)
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
      madata.kdtree_coords->sort_results = true;
      
   // Inside processing
   {
      sb_points(input_parameters, madata, 1);
#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Done shrinking interior balls, took " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif
   }

   // Outside processing
   {
      sb_points(input_parameters, madata, 0);
#ifdef VERBOSEPRINT
      t0.elapse();
      std::cout << "Done shrinking exterior balls, took " << t0.getTime()*1000.0 << " ms" << std::endl;
#endif
   }

   // Free memory
   delete madata.kdtree_coords; madata.kdtree_coords = NULL;
}

/*

void convertNPYtoXYZ(std::string inFile, std::string outFile)
{
   // Read in the data:
   cnpy::NpyArray coords_npy = cnpy::npy_load(inFile.c_str());
   float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

   unsigned int num_points = coords_npy.shape[0];
   unsigned int dim = coords_npy.shape[1];

   // Write this out to a pointcloudxyz file:
   outFile += ".xyz";
   std::ofstream out_pointcloudxyz(outFile);
   if (!out_pointcloudxyz)
      throw TCLAP::ArgParseException("invalid filepath", outFile);
   out_pointcloudxyz << "x y z\n";

   for (int i = 0; i < num_points; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         if (j > 0) out_pointcloudxyz << " ";
         out_pointcloudxyz << coords_carray[i*3 + j];
      }
      out_pointcloudxyz << "\n";
   }
   coords_npy.destruct();
}

*/

