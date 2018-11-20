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
#include "madata.h"
#include "compute_normals_processing.h"

using namespace masb;

int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Estimates normals using PCA, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' file; a Nx3 float array where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "output", "path to output directory. Estimated normals are written to the file 'normals.npy'.", false, "", "output dir", cmd);

        TCLAP::ValueArg<int> kArg("k","kneighbours","number of nearest neighbours to use for PCA",false,10,"int", cmd);

        TCLAP::SwitchArg reorder_kdtreeSwitch("N","no-kdtree-reorder","Don't reorder kd-tree points: slower computation but lower memory use", cmd, true);
        
        cmd.parse(argc,argv);
        
        normals_parameters input_parameters;
        input_parameters.k = kArg.getValue();

        input_parameters.kd_tree_reorder = reorder_kdtreeSwitch.getValue();

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();
        std::replace(output_path.begin(), output_path.end(), '\\', '/');


        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::replace(input_coords_path.begin(), input_coords_path.end(), '\\', '/');

        output_path += "/normals.npy";
        // check for proper in-output arguments
        {
            std::ifstream infile(input_coords_path.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ofstream outfile(output_path.c_str());    
            if(!outfile)
                throw TCLAP::ArgParseException("invalid filepath", output_path);
        }

        std::cout << "Parameters: k="<<input_parameters.k<<"\n";
        
        ma_data madata = {};

        cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
        float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

        madata.m = coords_npy.shape[0];
        unsigned int dim = coords_npy.shape[1];
        PointList coords(madata.m);
        for (unsigned int i=0; i<madata.m; i++) coords[i] = Point(&coords_carray[i*3]);
        coords_npy.destruct();

        VectorList normals(madata.m);
        madata.normals = &normals;
        madata.coords = &coords;

        // Perform the actual processing
        compute_normals(input_parameters, madata);

        // Output results
        Scalar* normals_carray = new Scalar[madata.m * 3];
        for (int i = 0; i < normals.size(); i++)
           for (int j = 0; j < 3; j++)
              normals_carray[i * 3 + j] = normals[i][j];

        const unsigned int c_size = (unsigned int) normals.size();
        const unsigned int shape[] = { c_size,3 };
        cnpy::npy_save(output_path.c_str(), normals_carray, shape, 2, "w");

        // Free memory
        delete[] normals_carray; normals_carray = NULL;
    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
