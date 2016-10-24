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

// tclap
#include <tclap/CmdLine.h>

// typedefs
#include "simplify_processing.h"
#include "io.h"



int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Feature-aware pointcloud simplification based on the Medial Axis Transform, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to input directory with inside it a 'coords.npy' and 'ma_*.npy' files. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory", false, "", "output dir", cmd);

        TCLAP::ValueArg<double> epsilonArg("e","epsilon","Control the degree of simplification, higher values mean more simplification. Typical values are in the range [0.01,0.6].",false,0.1,"double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c","cellsize","Cellsize used during grid-based lfs simplification (in units of your dataset). Large cellsize means faster processing, but potentially more noticable jumps in point density at cell boundaries.",false,1,"double", cmd);
        TCLAP::ValueArg<double> bisecArg("b","bisec","Bisector threshold used to clean the MAT points before LFS computation. With lower values more aggressive cleaning is performed which means more robustness to noise (in the MAT) but also less features will be detected. Typical range [0.1,10] (degrees).",false,1,"double", cmd);
        TCLAP::ValueArg<double> maxdensArg("m","max","Upper bound point density in pts/unit^2",false,1,"double", cmd);
        
        TCLAP::ValueArg<double> fake3dArg("f","fake3d","Use 2D grid instead of 3D grid, intended for 2.5D datasets (eg. buildings without points only on the roof and not on the walls). In addition this mode will try to detect elevation jumps in the dataset (eg. where there should be a wall) and still try to preserve points around those areas, the value for this parameter is the threshold elevation difference (in units of your dataset) within one gridcell that will be used for the elevation jump detection function.",false,0.5,"double", cmd);
        TCLAP::SwitchArg innerSwitch("i","inner","Compute LFS using only interior MAT points.", cmd, false);
        TCLAP::SwitchArg squaredSwitch("s","squared","Use squared LFS during simplification.", cmd, false);
        TCLAP::SwitchArg nolfsSwitch("d","no-lfs","Don't recompute lfs.'", cmd, false);
        
        TCLAP::ValueArg<std::string> outputXYZArg("a","xyz","output filtered points to plain .xyz text file",false,"lfs_simp.xyz","string", cmd);

        cmd.parse(argc,argv);
        
        simplify_parameters input_parameters;

        input_parameters.epsilon = epsilonArg.getValue();
        input_parameters.cellsize = cellsizeArg.getValue();
        input_parameters.bisec_threshold = (bisecArg.getValue() / 180.0) * M_PI;
        
        input_parameters.compute_lfs = !nolfsSwitch.getValue();
        input_parameters.elevation_threshold = fake3dArg.getValue();
        input_parameters.maximum_density = maxdensArg.getValue();
        input_parameters.true_z_dim = true;
        input_parameters.only_inner = innerSwitch.getValue();
        input_parameters.squared = squaredSwitch.getValue();
        if( fake3dArg.isSet() )
           input_parameters.true_z_dim = false;

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();
        std::replace(output_path.begin(), output_path.end(), '\\', '/');


        ma_data madata = {};
        io_parameters input_params = {};
        input_params.coords = true;
        input_params.ma_coords = true;
        input_params.ma_qidx = true;
        if(!input_parameters.compute_lfs){
           input_params.lfs = true;
        }
        else
        {
           madata.lfs.resize(madata.m);
        }

        npy2madata(inputArg.getValue(), madata, input_params);

        madata.mask.resize(madata.m);

	    {
          // Perform the actual processing
          simplify_lfs(input_parameters, madata);
          
          // count number of remaining points
          unsigned int cnt;
          for( int i=0; i<madata.m; i++ )
                if( madata.mask[i] ) cnt++;
          std::cout << cnt << " out of " << madata.m << " points remaining [" << int(100*float(cnt)/madata.m) << "%]" << std::endl;

          // Output results
          io_parameters output_params = {};
          output_params.lfs = true;
          output_params.mask = true;
          madata2npy(output_path, madata, output_params);
        }

        if( outputXYZArg.isSet() ){
            /*
            std::string outFile_bounds = outputXYZArg.getValue();
            outFile_bounds.append(".bounds");
            std::replace(outFile_bounds.begin(), outFile_bounds.end(), '\\', '/');
            
            std::ofstream ofs_bounds(outFile_bounds.c_str());
            ofs_bounds << madata.bbox.min[0] << std::endl << madata.bbox.max[0] << 
             std::endl << madata.bbox.min[1] << std::endl << madata.bbox.max[1] << 
             std::endl << madata.bbox.min[2] << std::endl << madata.bbox.max[2] << std::endl;
            
            ofs_bounds.close();
            */

            std::string outFile_xyz = outputXYZArg.getValue();
            std::replace(outFile_xyz.begin(), outFile_xyz.end(), '\\', '/');
            std::ofstream ofs(outFile_xyz.c_str());
            //ofs <<std::setprecision(2)<<std::fixed;
            
            // many pointcloud xyz readers prefer a "header" line.
            ofs << "x y z" << std::endl;

            for( int i=0; i<madata.m; i++ ) {
                if( madata.mask[i] ){
                    ofs << (*madata.coords)[i].x;
                    ofs << " " << (*madata.coords)[i].y;
                    ofs << " " << (*madata.coords)[i].z;
                ofs << std::endl;
                }
            }
            
            ofs.close();
        }
	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
