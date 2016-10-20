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

#include <iostream>
#include <fstream>
#include <string>

#include <tclap/CmdLine.h>

#include "compute_normals_processing.h"
#include "io.h"
#include "madata.h"
#include "types.h"

int main(int argc, char **argv) {
   // parse command line arguments
   try {
      TCLAP::CmdLine cmd("Estimates normals using PCA, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

      TCLAP::UnlabeledValueArg<std::string> inputArg("input", "path to directory with inside it a 'coords.npy' file; a Nx3 float array where N is the number of input points.", true, "", "input dir", cmd);
      TCLAP::UnlabeledValueArg<std::string> outputArg("output", "path to output directory. Estimated normals are written to the file 'normals.npy'.", false, "", "output dir", cmd);

      TCLAP::ValueArg<int> kArg("k", "kneighbours", "number of nearest neighbours to use for PCA", false, 10, "int", cmd);

      TCLAP::SwitchArg reorder_kdtreeSwitch("N", "no-kdtree-reorder", "Don't reorder kd-tree points: slower computation but lower memory use", cmd, true);

      cmd.parse(argc, argv);

      normals_parameters normal_params;
      normal_params.k = kArg.getValue();
      normal_params.kd_tree_reorder = reorder_kdtreeSwitch.getValue();

      std::string npy_path = inputArg.getValue();
      std::string npy_path_coords = npy_path + "/coords.npy";

      if (outputArg.isSet()) {
         npy_path = outputArg.getValue();
      }

      std::cout << "Parameters: k=" << normal_params.k << std::endl;

      io_parameters io_params = {};
      io_params.coords = true;
      ma_data madata = {};

      PointCloud::Ptr coords(new PointCloud);
      madata.coords = coords;
      NormalCloud::Ptr normals(new NormalCloud);
      madata.normals = normals;

      npy2madata(npy_path, madata, io_params);

      std::cout << "Point count: " << madata.m << std::endl;

      // Perform the actual processing
      compute_normals(normal_params, madata);

      io_params.coords = false;
      io_params.normals = true;
      madata2npy(npy_path, madata, io_params);
   }
   catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

   return 0;
}
