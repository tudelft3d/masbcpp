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
#include "simplify_processing.h"



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
        input_parameters.dimension = 3;
        input_parameters.only_inner = innerSwitch.getValue();
        input_parameters.squared = squaredSwitch.getValue();
        if( fake3dArg.isSet() )
           input_parameters.dimension = 2;

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();
        std::replace(output_path.begin(), output_path.end(), '\\', '/');

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_path = inputArg.getValue();
        std::replace(input_path.begin(), input_path.end(), '\\', '/');
        std::string input_coords_path = input_path+"/coords.npy";
        std::string input_path_ma_coords_in = input_path+"/ma_coords_in.npy";
        std::string input_path_ma_coords_out = input_path+"/ma_coords_out.npy";
        std::string input_path_ma_qidx_in = input_path+"/ma_qidx_in.npy";
        std::string input_path_ma_qidx_out = input_path+"/ma_qidx_out.npy";
        std::string input_path_lfs = input_path+"/lfs.npy";
        std::string output_lfs = output_path+"/lfs.npy";
        std::string output_filtermask = output_path+"/decimate_lfs.npy";
        {
            std::ifstream infile(input_coords_path.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ifstream infile(input_path_ma_coords_in.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ofstream outfile(output_filtermask.c_str());    
            if(!outfile)
                throw TCLAP::ArgParseException("invalid filepath", output_path);
        }

        ma_data madata = {};
	    
	    cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
	    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

	    madata.m = coords_npy.shape[0];
	    unsigned int dim = coords_npy.shape[1];
	    PointList coords(madata.m);
        madata.bbox = Box(Point(&coords_carray[0]), Point(&coords_carray[0]));
	    for ( int i=0; i<madata.m; i++){ 
            coords[i] = Point(&coords_carray[i*3]);
            madata.bbox.addPoint(coords[i]);
        }
	    coords_npy.destruct();
        // std::cout << "bbox: " << madata.bbox.min[0] << " " << madata.bbox.min[1] << " " << madata.bbox.min[2] << " " << madata.bbox.max[0] << " " << madata.bbox.max[1] << " " << madata.bbox.max[2] << std::endl;

        PointList ma_coords(2*madata.m);

        cnpy::NpyArray ma_coords_in_npy = cnpy::npy_load( input_path_ma_coords_in.c_str() );
        float* ma_coords_in_carray = reinterpret_cast<float*>(ma_coords_in_npy.data);
        for ( int i=0; i<madata.m; i++) ma_coords[i] = Point(&ma_coords_in_carray[i*3]);
        ma_coords_in_npy.destruct();    
        
        cnpy::NpyArray ma_coords_out_npy = cnpy::npy_load( input_path_ma_coords_out.c_str() );
        float* ma_coords_out_carray = reinterpret_cast<float*>(ma_coords_out_npy.data);
        for ( int i=0; i<madata.m; i++) ma_coords[i+madata.m] = Point(&ma_coords_out_carray[i*3]);
        ma_coords_out_npy.destruct();

        
        madata.ma_qidx = new int[2*madata.m];
        cnpy::NpyArray ma_qidx_in_npy = cnpy::npy_load( input_path_ma_qidx_in.c_str() );
        int* ma_qidx_in = reinterpret_cast<int*>(ma_qidx_in_npy.data);
        for ( int i=0; i<madata.m; i++) madata.ma_qidx[i] = ma_qidx_in[i];
        ma_qidx_in_npy.destruct();

        cnpy::NpyArray ma_qidx_out_npy = cnpy::npy_load( input_path_ma_qidx_out.c_str() );
        int* ma_qidx_out = reinterpret_cast<int*>(ma_qidx_out_npy.data);
        for ( int i=0; i<madata.m; i++) madata.ma_qidx[i+madata.m] = ma_qidx_out[i];
        ma_qidx_out_npy.destruct();

        if(input_parameters.compute_lfs){
            madata.lfs = new float[madata.m];
        } else {
            cnpy::NpyArray ma_lfs_npy = cnpy::npy_load( input_path_lfs.c_str() );
            madata.lfs = reinterpret_cast<float*>(ma_lfs_npy.data);
        }
	    
        
        madata.coords = &coords; // don't own this memory
        // madata.normals = &normals;
        madata.ma_coords = &ma_coords; // don't own this memory

        madata.mask = new bool[madata.m];

	    {
          // Perform the actual processing
          simplify_lfs(input_parameters, madata);
          
          // count number of remaining points
          unsigned int cnt;
          for( int i=0; i<madata.m; i++ )
                if( madata.mask[i] ) cnt++;
          std::cout << cnt << " out of " << madata.m << " points remaining [" << int(100*float(cnt)/madata.m) << "%]" << std::endl;

          // Output results
          const unsigned int c_size = madata.m;
          const unsigned int shape[] = { c_size };
          cnpy::npy_save(output_filtermask.c_str(), madata.mask, shape, 1, "w");
          cnpy::npy_save(output_lfs.c_str(), madata.lfs, shape, 1, "w");
        }

        if( outputXYZArg.isSet() ){
            std::string outFile_bounds = outputXYZArg.getValue();
            outFile_bounds.append(".bounds");
            std::replace(outFile_bounds.begin(), outFile_bounds.end(), '\\', '/');
            
            std::ofstream ofs_bounds(outFile_bounds.c_str());
            ofs_bounds << madata.bbox.min[0] << std::endl << madata.bbox.max[0] << 
             std::endl << madata.bbox.min[1] << std::endl << madata.bbox.max[1] << 
             std::endl << madata.bbox.min[2] << std::endl << madata.bbox.max[2] << std::endl;
            
            ofs_bounds.close();
            
            outFile_bounds = outputXYZArg.getValue();
            std::replace(outFile_bounds.begin(), outFile_bounds.end(), '\\', '/');
            std::ofstream ofs(outFile_bounds.c_str());
            ofs <<std::setprecision(2)<<std::fixed;
            
            // many pointcloud xyz readers prefer a "header" line.
            ofs << "x y z" << std::endl;

            for( int i=0; i<madata.m; i++ ) {
                if( madata.mask[i] ){
                    ofs << (*madata.coords)[i][0];
                    ofs << " " << (*madata.coords)[i][1];
                    ofs << " " << (*madata.coords)[i][2];
                ofs << std::endl;
                }
            }
            
            ofs.close();
        }

        // Free memory
        delete[] madata.mask; madata.mask = NULL;
        delete[] madata.lfs; madata.lfs = NULL;
        delete[] madata.ma_qidx; madata.ma_qidx = NULL;

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
