// Copyright (c) 2015
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// This file is part of masbcpp.
//
// masbcpp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// masbcpp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with masbcpp.  If not, see <http://www.gnu.org/licenses/>.

// #define VERBOSEPRINT 1;
// #define WITH_OPENMP 1;

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <cmath>

// OpenMP
#ifdef WITH_OPENMP
    #include <omp.h>
#endif

// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>
#ifndef __MINGW32__
#include <vrui/Misc/Timer.h>
#endif
// kdtree2
#include <kdtree2/kdtree2.hpp>
// cnpy
#include <cnpy/cnpy.h>
// tclap
#include <tclap/CmdLine.h>

// typedefs
#include "types.h"

// globals
const Point nanPoint( std::numeric_limits<Scalar>::quiet_NaN() );


compute_lfs(ma_data &madata)
{
    #ifndef __MINGW32__
    Misc::Timer t0;
    #endif
    kdtree2::KDTree* kd_tree;
    kd_tree = new kdtree2::KDTree(madata.ma_coords_in,True);
    kd_tree->sort_results = true;
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Constructed kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif

    kdtree2::KDTreeResultVector result;
    for( unsigned int i=0; i<madata.m; i++ ){
        kd_tree->n_nearest(madata.points[i], 2, result);
        madata.lfs[i] = sqrt(result[1].dis);
    }
}


void sb_points(PointList &points, VectorList &normals, kdtree2::KDTree* kd_tree, PointList &ma_coords, int* ma_qidx, bool inner=1)
{
    Point p;
    Vector n;

    #pragma omp parallel for private(p, n)
    for( unsigned int i=0; i<points.size(); i++ )
    {
        p = points[i];
        if( inner )
            n = normals[i];
        else
            n = -normals[i];
        ma_result r = sb_point(p, n, kd_tree);
        ma_coords[i] = r.c;
        ma_qidx[i] = r.qidx;
    }
    // return ma_coords;
}


int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory. Estimated MAT points are written to the files 'ma_coords_in.npy' and 'ma_coords_out.npy', for interior and exterior MAT respectively.", true, "", "output dir", cmd);

        TCLAP::ValueArg<double> epsiolonArg("e","preserve","denoise preserve threshold",false,20,"double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c","planar","denoise planar threshold",false,32,"double", cmd);
        
        TCLAP::SwitchArg towdimSwitch("t","twodim","use 2D grid instead of 3D grid", cmd, false);

        cmd.parse(argc,argv);
        
        epsilon = epsilonArg.getValue();
        cellsize = cellsizeArg.getValue();
        
        bool twodim = towdimSwitch.getValue();

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_path_ma_in = inputArg.getValue()+"/ma_coords_in.npy";
        std::string input_path_ma_out = inputArg.getValue()+"/ma_coords_out.npy";
        std::string output_filtermask = outputArg.getValue()+"/decimate_lfs.npy";
        {
            std::ifstream infile(input_coords_path.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ifstream infile(input_path_ma_in.c_str());
            if(!infile)
                throw TCLAP::ArgParseException("invalid filepath", inputArg.getValue());
        }
        {
            std::ofstream outfile(output_filtermask.c_str());    
            if(!outfile)
                throw TCLAP::ArgParseException("invalid filepath", outputArg.getValue());
        }
	    
	    cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
	    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

	    unsigned int num_points = coords_npy.shape[0];
	    unsigned int dim = coords_npy.shape[1];
	    PointList coords(num_points);
	    for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
	    coords_npy.destruct();

        cnpy::NpyArray ma_coords_in_npy = cnpy::npy_load( input_ma_coords_in_path.c_str() );
        float* ma_coords_in_carray = reinterpret_cast<float*>(ma_coords_in_npy.data);
        VectorList ma_coords_in(ma_coords_in_npy.shape[0]);
        for ( int i=0; i<num_points; i++) ma_coords_in[i] = Vector(&ma_coords_in_carray[i*3]);
        ma_coords_in_npy.destruct();
        
	    cnpy::NpyArray ma_coords_out_npy = cnpy::npy_load( input_ma_coords_out_path.c_str() );
	    float* ma_coords_out_carray = reinterpret_cast<float*>(ma_coords_out_npy.data);
	    VectorList ma_coords_in(ma_coords_out_npy.shape[0]);
	    for ( int i=0; i<num_points; i++) ma_coords_in[i] = Vector(&ma_coords_out_carray[i*3]);
	    ma_coords_out_npy.destruct();
	    
        ma_data madata = {};
        madata.m = num_points;
        madata.coords = &coords;
        madata.normals = &normals;
        madata.ma_coords_in = &ma_coords_in;
        madata.ma_coords_out = &ma_coords_out;

        madata.mask = new bool[num_points];
        madata.lfs = new float[num_points];

	    // omp_set_num_threads(4);

	    {
            bool* madata.mask

            // compute lfs
            compute_lfs(ma_data &madata)
            // create filter ...
	    
	        const unsigned int c_size = num_points
            const unsigned int shape_[] = {c_size};
	        cnpy::npy_save(output_filtermask.c_str(), madata.mask, shape_, 1, "w");
	    }

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
