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

#include "compute_ma_processing.h"
#include "io.h"
#include "madata.h"
#include "types.h"

int main(int argc, char **argv) {
   // parse command line arguments
   try {
      TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

      TCLAP::UnlabeledValueArg<std::string> inputArg("input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
      TCLAP::UnlabeledValueArg<std::string> outputArg("output", "path to output directory", false, "", "output dir", cmd);

      TCLAP::ValueArg<double> denoise_preserveArg("d", "preserve", "denoise preserve threshold", false, 20, "double", cmd);
      TCLAP::ValueArg<double> denoise_planarArg("p", "planar", "denoise planar threshold", false, 32, "double", cmd);
      TCLAP::ValueArg<double> initial_radiusArg("r", "radius", "initial ball radius", false, 200, "double", cmd);

      TCLAP::SwitchArg nan_for_initrSwitch("a", "nan", "write nan for points with radius equal to initial radius", cmd, false);

      cmd.parse(argc, argv);

      ma_parameters input_parameters;

      input_parameters.initial_radius = float(initial_radiusArg.getValue());
      input_parameters.denoise_preserve = (M_PI / 180.0) * denoise_preserveArg.getValue();
      input_parameters.denoise_planar = (M_PI / 180.0) * denoise_planarArg.getValue();
      input_parameters.nan_for_initr = nan_for_initrSwitch.getValue();

      std::string output_path = inputArg.getValue();
      if (outputArg.isSet())
         output_path = outputArg.getValue();
      std::replace(output_path.begin(), output_path.end(), '\\', '/');

      // check for proper in-output arguments and set in and output filepath strings
      std::string input_coords_path = inputArg.getValue() + "/coords.npy";
      std::replace(input_coords_path.begin(), input_coords_path.end(), '\\', '/');
      std::string input_normals_path = inputArg.getValue() + "/normals.npy";
      std::replace(input_normals_path.begin(), input_normals_path.end(), '\\', '/');
      std::string output_path_ma_in = output_path + "/ma_coords_in.npy";
      std::string output_path_ma_out = output_path + "/ma_coords_out.npy";
      std::string output_path_ma_q_in = output_path + "/ma_qidx_in.npy";
      std::string output_path_ma_q_out = output_path + "/ma_qidx_out.npy";

      std::string output_path_metadata = output_path + "/compute_ma";
      std::replace(output_path_metadata.begin(), output_path_metadata.end(), '\\', '/');

      {
         std::ifstream infile(input_coords_path.c_str());
         if (!infile)
            throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
      }
      {
         std::ifstream infile(input_normals_path.c_str());
         if (!infile)
            throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
      }
      {
         std::ofstream outfile(output_path_ma_in.c_str());
         if (!outfile)
            throw TCLAP::ArgParseException("invalid filepath", output_path);
      }

      std::cout << "Parameters: denoise_preserve=" << denoise_preserveArg.getValue() << ", denoise_planar=" << denoise_planarArg.getValue() << ", initial_radius=" << input_parameters.initial_radius << "\n";

      io_parameters io_params = {};
      io_params.coords = true;
      io_params.normals = true;
      ma_data madata = {};

      madata.coords = PointCloud::Ptr(new PointCloud);
      madata.normals = NormalCloud::Ptr(new NormalCloud);
      madata.ma_coords = PointCloud::Ptr(new PointCloud);

      npy2madata(inputArg.getValue(), madata, io_params);

      madata.ma_qidx.resize(2 * madata.m);

      // Perform the actual processing
      compute_masb_points(input_parameters, madata);

      io_params.coords = false;
      io_params.normals = false;
      io_params.ma_coords = true;
      io_params.ma_qidx = true;
      madata2npy(output_path, madata, io_params);
   }
   catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

   return 0;
}
