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

// cnpy
#include <cnpy/cnpy.h>
// tclap
#include <tclap/CmdLine.h>

// typedefs
#include "types.h"
#include "compute_ma_processing.h"

int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "output", "path to output directory", false, "", "output dir", cmd);

        TCLAP::ValueArg<double> denoise_preserveArg("d","preserve","denoise preserve threshold",false,20,"double", cmd);
        TCLAP::ValueArg<double> denoise_planarArg("p","planar","denoise planar threshold",false,32,"double", cmd);
        TCLAP::ValueArg<double> initial_radiusArg("r","radius","initial ball radius",false,200,"double", cmd);
        
        TCLAP::SwitchArg nan_for_initrSwitch("a","nan","write nan for points with radius equal to initial radius", cmd, false);
        TCLAP::SwitchArg reorder_kdtreeSwitch("N","no-kdtree-reorder","Don't reorder kd-tree points: slower computation but lower memory use", cmd, true);

        cmd.parse(argc,argv);
        
        ma_parameters input_parameters;

        input_parameters.initial_radius = initial_radiusArg.getValue();
        input_parameters.denoise_preserve = (M_PI/180) * denoise_preserveArg.getValue();
        input_parameters.denoise_planar = (M_PI/180) * denoise_planarArg.getValue();
        
        input_parameters.nan_for_initr = nan_for_initrSwitch.getValue();
        input_parameters.kd_tree_reorder = reorder_kdtreeSwitch.getValue();

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();
        std::replace(output_path.begin(), output_path.end(), '\\', '/');

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::replace(input_coords_path.begin(), input_coords_path.end(), '\\', '/');
        std::string input_normals_path = inputArg.getValue()+"/normals.npy";
        std::replace(input_normals_path.begin(), input_normals_path.end(), '\\', '/');
        std::string output_path_ma_in = output_path+"/ma_coords_in.npy";
        std::string output_path_ma_out = output_path+"/ma_coords_out.npy";
        std::string output_path_ma_q_in = output_path+"/ma_qidx_in.npy";
        std::string output_path_ma_q_out = output_path+"/ma_qidx_out.npy";
        {
            std::ifstream infile(input_coords_path.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ifstream infile(input_normals_path.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ofstream outfile(output_path_ma_in.c_str());    
            if(!outfile)
                throw TCLAP::ArgParseException("invalid filepath", output_path);
        }

	   	std::cout << "Parameters: denoise_preserve="<<denoise_preserveArg.getValue()<<", denoise_planar="<<denoise_planarArg.getValue()<<", initial_radius="<<input_parameters.initial_radius<<"\n";
	    
	    cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
	    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

	    unsigned int num_points = coords_npy.shape[0];
	    unsigned int dim = coords_npy.shape[1];
	    PointList coords(num_points);
	    for ( unsigned int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
	    coords_npy.destruct();

	    cnpy::NpyArray normals_npy = cnpy::npy_load( input_normals_path.c_str() );
	    float* normals_carray = reinterpret_cast<float*>(normals_npy.data);
	    VectorList normals(normals_npy.shape[0]);
	    for ( unsigned int i=0; i<num_points; i++) normals[i] = Vector(&normals_carray[i*3]);
	    normals_npy.destruct();

       // Storage space for our results:
       PointList ma_coords_in(coords.size());
       int* ma_qidx_in = new int[num_points];
       PointList ma_coords_out(coords.size());
       int* ma_qidx_out = new int[num_points];

       // Perform the actual processing
       compute_masb_points(input_parameters, coords, normals, ma_coords_in, ma_qidx_in, ma_coords_out, ma_qidx_out);

       // Write out the results for the inside
       Scalar* ma_coords_in_carray = new Scalar[num_points * 3];
       for (int i = 0; i < ma_coords_in.size(); i++)
          for (int j = 0; j < 3; j++)
             ma_coords_in_carray[i * 3 + j] = ma_coords_in[i][j];

       const unsigned int c_size_in = ma_coords_in.size();
       const unsigned int shape_in[] = { c_size_in,3 };
       cnpy::npy_save(output_path_ma_in.c_str(), ma_coords_in_carray, shape_in, 2, "w");
       const unsigned int shape_in_[] = { c_size_in };
       cnpy::npy_save(output_path_ma_q_in.c_str(), ma_qidx_in, shape_in_, 1, "w");

       // Write out the results for the outside
       Scalar* ma_coords_out_carray = new Scalar[num_points * 3];
       for (int i = 0; i < ma_coords_out.size(); i++)
          for (int j = 0; j < 3; j++)
             ma_coords_out_carray[i * 3 + j] = ma_coords_out[i][j];

       const unsigned int c_size_out = ma_coords_out.size();
       const unsigned int shape_out[] = { c_size_out,3 };
       cnpy::npy_save(output_path_ma_out.c_str(), ma_coords_out_carray, shape_in, 2, "w");
       const unsigned int shape_out_[] = { c_size_out };
       cnpy::npy_save(output_path_ma_q_out.c_str(), ma_qidx_out, shape_out_, 1, "w");

       // Free up memory
       delete[] ma_qidx_in; ma_qidx_in = NULL;
       delete[] ma_qidx_out; ma_qidx_out = NULL;

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
