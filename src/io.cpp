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

#include "io.h"

#include <iostream>
#include <fstream>
#include <string>

#include <cnpy/cnpy.h>

#include "madata.h"
#include "types.h"

inline cnpy::NpyArray read_npyarray(std::string input_file_path) {
   // windows fix
   std::replace(input_file_path.begin(), input_file_path.end(), '\\', '/');
   // check if file exists
   {
      std::ifstream infile(input_file_path.c_str());
      if (!infile) {
         std::cerr << "Invalid file path " << input_file_path << std::endl;
         exit(1);
      }
   }
   // std::cout << "Reading array from " << input_file_path <<std::endl;
   cnpy::NpyArray npy_array = cnpy::npy_load(input_file_path.c_str());
   return npy_array;
}

void npy2madata(std::string input_dir_path, ma_data &madata, io_parameters &params) {
   if (params.coords) {
      std::cout << "Reading coords array..." << std::endl;

      cnpy::NpyArray npy_array = read_npyarray(input_dir_path + "/coords.npy");
      float* coords_carray = reinterpret_cast<float*>(npy_array.data);

      madata.coords.reset(new PointCloud);
      madata.coords->reserve(npy_array.shape[0]);

      for (size_t i = 0; i < npy_array.shape[0]; i++)
         madata.coords->push_back(Point(
            coords_carray[i * 3 + 0],
            coords_carray[i * 3 + 1],
            coords_carray[i * 3 + 2]
         ));
      npy_array.destruct();
   }

   if (params.normals) {
      std::cout << "Reading normals array..." << std::endl;

      cnpy::NpyArray npy_array = read_npyarray(input_dir_path + "/normals.npy");
      float* normals_carray = reinterpret_cast<float*>(npy_array.data);

      if (npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and normals" << std::endl;
         exit(1);
      }

      madata.normals.reset(new NormalCloud);
      madata.normals->reserve(madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.normals->push_back(Normal(
            normals_carray[i * 3 + 0],
            normals_carray[i * 3 + 1],
            normals_carray[i * 3 + 2]
         ));
      npy_array.destruct();
   }

   if (params.ma_coords) {
      std::cout << "Reading ma coords arrays..." << std::endl;

      cnpy::NpyArray in_npy_array = read_npyarray(input_dir_path + "/ma_coords_in.npy");
      float* in_ma_coords_carray = reinterpret_cast<float*>(in_npy_array.data);

      if (in_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and inner ma coords" << std::endl;
         exit(1);
      }

      cnpy::NpyArray out_npy_array = read_npyarray(input_dir_path + "/ma_coords_out.npy");
      float* out_ma_coords_carray = reinterpret_cast<float*>(out_npy_array.data);

      if (out_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and outer ma coords" << std::endl;
         exit(1);
      }

      madata.ma_coords.reset(new PointCloud);
      madata.ma_coords->reserve(2 * madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_coords->push_back(Point(
            in_ma_coords_carray[i * 3 + 0],
            in_ma_coords_carray[i * 3 + 1],
            in_ma_coords_carray[i * 3 + 2]
         ));
      in_npy_array.destruct();

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_coords->push_back(Point(
            out_ma_coords_carray[i * 3 + 0],
            out_ma_coords_carray[i * 3 + 1],
            out_ma_coords_carray[i * 3 + 2]
         ));
      out_npy_array.destruct();
   }

   if (params.ma_qidx) {
      std::cout << "Reading q index arrays..." << std::endl;

      cnpy::NpyArray in_npy_array = read_npyarray(input_dir_path + "/ma_qidx_in.npy");
      int* in_qidx_carray = reinterpret_cast<int*>(in_npy_array.data);

      if (in_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and inner q indices" << std::endl;
         exit(1);
      }

      cnpy::NpyArray out_npy_array = read_npyarray(input_dir_path + "/ma_qidx_out.npy");
      int* out_qidx_carray = reinterpret_cast<int*>(out_npy_array.data);

      if (out_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and outer q indices" << std::endl;
         exit(1);
      }

      madata.ma_qidx.reserve(2 * madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_qidx.push_back(in_qidx_carray[i]);
      in_npy_array.destruct();

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_qidx.push_back(out_qidx_carray[i]);
      out_npy_array.destruct();
   }

   if (params.lfs) {
      std::cout << "Reading lfs array..." << std::endl;

      cnpy::NpyArray npy_array = read_npyarray(input_dir_path + "/lfs.npy");
      float* lfs_carray = reinterpret_cast<float*>(npy_array.data);

      if (npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and lfs" << std::endl;
         exit(1);
      }

      madata.lfs.reserve(madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.lfs.push_back(lfs_carray[i]);
      npy_array.destruct();
   }
}

void madata2npy(std::string npy_path, ma_data &madata, io_parameters &params) {
   if (params.coords) {
      std::cout << "Writing coords array..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };
      float* coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         coords_carray[i * 3 + 0] = madata.coords->at(i).x;
         coords_carray[i * 3 + 1] = madata.coords->at(i).y;
         coords_carray[i * 3 + 2] = madata.coords->at(i).z;
      }
      cnpy::npy_save(npy_path + "/coords.npy", coords_carray, shape, 2, "w");
      delete[] coords_carray; coords_carray = nullptr;
   }

   if (params.normals) {
      std::cout << "Writing normals array..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };
      float* normals_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         normals_carray[i * 3 + 0] = madata.normals->at(i).normal_x;
         normals_carray[i * 3 + 1] = madata.normals->at(i).normal_y;
         normals_carray[i * 3 + 2] = madata.normals->at(i).normal_z;
      }
      cnpy::npy_save(npy_path + "/normals.npy", normals_carray, shape, 2, "w");
      delete[] normals_carray; normals_carray = nullptr;
   }

   if (params.ma_coords) {
      std::cout << "Writing ma coords arrays..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };

      float* in_ma_coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         in_ma_coords_carray[i * 3 + 0] = madata.ma_coords->at(i).x;
         in_ma_coords_carray[i * 3 + 1] = madata.ma_coords->at(i).y;
         in_ma_coords_carray[i * 3 + 2] = madata.ma_coords->at(i).z;
      }
      cnpy::npy_save(npy_path + "/ma_coords_in.npy", in_ma_coords_carray, shape, 2, "w");
      delete[] in_ma_coords_carray; in_ma_coords_carray = nullptr;

      float* out_ma_coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         out_ma_coords_carray[i * 3 + 0] = madata.ma_coords->at(i + madata.coords->size()).x;
         out_ma_coords_carray[i * 3 + 1] = madata.ma_coords->at(i + madata.coords->size()).y;
         out_ma_coords_carray[i * 3 + 2] = madata.ma_coords->at(i + madata.coords->size()).z;
      }
      cnpy::npy_save(npy_path + "/ma_coords_out.npy", out_ma_coords_carray, shape, 2, "w");
      delete[] out_ma_coords_carray; out_ma_coords_carray = nullptr;
   }

   if (params.ma_qidx) {
      std::cout << "Writing q index arrays..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()) };

      cnpy::npy_save(npy_path + "/ma_qidx_in.npy", &madata.ma_qidx[0], shape, 1, "w");
      cnpy::npy_save(npy_path + "/ma_qidx_out.npy", &madata.ma_qidx[madata.coords->size()], shape, 1, "w");
   }

   if (params.lfs) {
      std::cout << "Writing lfs array..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()) };
      cnpy::npy_save(npy_path + "/lsf.npy", &madata.lfs[0], shape, 1, "w");
   }

   if (params.mask) {
      std::cout << "Writing mask array..." << std::endl;

      const unsigned int shape[] = { static_cast<unsigned int>(2*madata.coords->size()) };
      float* out_mask_carray = new float[2*madata.coords->size()];
      for (size_t i = 0; i < 2*madata.coords->size(); i++) {
         out_mask_carray[i] = madata.mask[i];
      }
      cnpy::npy_save(npy_path + "/decimate_lfs.npy", out_mask_carray, shape, 1, "w");
      delete[] out_mask_carray; out_mask_carray = nullptr;
   }
}
