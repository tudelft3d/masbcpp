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


void compute_lfs(ma_data &madata, float bisec_threshold)
{
    #ifndef __MINGW32__
    Misc::Timer t0;
    #endif

    // compute bisectors
    VectorList ma_bisec_in(madata.m);
    madata.ma_bisec_in = &ma_bisec_in;
    VectorList ma_bisec_out(madata.m);
    madata.ma_bisec_out = &ma_bisec_out;
    for( int i=0; i<madata.m; i++ ){

        Vector f1_in = (*madata.coords)[i] - (*madata.ma_coords_in)[i];
        Vector f2_in = (*madata.coords)[ madata.ma_qidx_in[i] ] - (*madata.ma_coords_in)[i];
        // madata.lfs[i] = f1_in[0];
        // std::cout << "ma_coords_in: " << (*madata.ma_coords_in)[i][0] << " " << (*madata.ma_coords_in)[i][1] << " " << (*madata.ma_coords_in)[i][2] << std::endl;
        // std::cout << "coords: " << (*madata.coords)[i][0] << " " << (*madata.coords)[i][1] << " " << (*madata.coords)[i][2] << std::endl;
        // std::cout << "f2_in: " << f2_in[0] << " " << f2_in[1] << " " << f2_in[2] << std::endl;
        // std::cout << "ma_qidx_in: " << madata.ma_qidx_in[i] << std::endl;

//        f1_in = f1_in.normalize();
//        f2_in = f2_in.normalize();

        (*madata.ma_bisec_in)[i] = (f1_in+f2_in).normalize();
        // std::cout << "ma_bisec_in: " << (*madata.ma_bisec_in)[i][0] << " " << (*madata.ma_bisec_in)[i][1] << " " << (*madata.ma_bisec_in)[i][2] << std::endl;

        Vector f1_out = (*madata.coords)[i] - (*madata.ma_coords_out)[i];
        Vector f2_out = (*madata.coords)[ madata.ma_qidx_out[i] ] - (*madata.ma_coords_out)[i];

//        f1_out = f1_out.normalize();
//        f2_out = f2_out.normalize();

        (*madata.ma_bisec_out)[i] = (f1_out+f2_out).normalize();

        // madata.ma_theta_in = np.arccos(np.sum(f1_in*f2_in,axis=1))
    }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Computed bisectors in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif

    int k = 2, count=0;
    {
        kdtree2::KDTree kd_tree(*madata.ma_coords_in,true);
        kd_tree.sort_results = true;
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Constructed kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
        
        kdtree2::KDTreeResultVector result;
        for( unsigned int i=0; i<madata.m; i++ ){
            kd_tree.n_nearest((*madata.ma_coords_in)[i], k, result);
            // compute bisector and filter .. rebuild kdtree .. compute lfs
            madata.mask[i] = false;
            // for( int j=1; j<k; j++ ){
                float bisec_angle = (*madata.ma_bisec_in)[result[1].idx] * (*madata.ma_bisec_in)[i];
                madata.lfs[i] = acos(bisec_angle);
                if( bisec_angle > bisec_threshold )
                    madata.mask[i] = true;
            // }
            if (madata.mask[i])
                count++;

        }
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Cleaned MA points in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
    }
    // mask and copy pointlist ma_coords
    PointList ma_coords_in_masked(count);
    int j=0;
    for( int i=0; i<madata.m; i++ ){
        if( madata.mask[i] )
            ma_coords_in_masked[j++] = (*madata.ma_coords_in)[i];
    }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Copied cleaned MA points in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif

    {
        // rebuild kd-tree
        kdtree2::KDTree kd_tree(ma_coords_in_masked,true);
        // kd_tree = new kdtree2::KDTree;
        kd_tree.sort_results = true;
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Constructed cleaned kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif

        kdtree2::KDTreeResultVector result;
        for( unsigned int i=0; i<madata.m; i++ ){
            kd_tree.n_nearest((*madata.coords)[i], k, result);
            // madata.lfs[i] = sqrt(result[1].dis);
        }
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Computed LFS "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
    }
}


int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory. Estimated MAT points are written to the files 'ma_coords_in.npy' and 'ma_coords_out.npy', for interior and exterior MAT respectively.", true, "", "output dir", cmd);

        TCLAP::ValueArg<double> epsilonArg("e","preserve","denoise preserve threshold",false,20,"double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c","planar","denoise planar threshold",false,32,"double", cmd);
        
        TCLAP::SwitchArg towdimSwitch("t","twodim","use 2D grid instead of 3D grid", cmd, false);

        cmd.parse(argc,argv);
        
        float epsilon = epsilonArg.getValue();
        float cellsize = cellsizeArg.getValue();
        
        bool twodim = towdimSwitch.getValue();

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_path_ma_coords_in = inputArg.getValue()+"/ma_coords_in.npy";
        std::string input_path_ma_coords_out = inputArg.getValue()+"/ma_coords_out.npy";
        std::string input_path_ma_qidx_in = inputArg.getValue()+"/ma_qidx_out.npy";
        std::string input_path_ma_qidx_out = inputArg.getValue()+"/ma_qidx_out.npy";
        std::string output_lfs = outputArg.getValue()+"/lfs.npy";
        std::string output_filtermask = outputArg.getValue()+"/decimate_lfs.npy";
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
                throw TCLAP::ArgParseException("invalid filepath", outputArg.getValue());
        }

        ma_data madata = {};
	    
	    cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
	    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

	    unsigned int num_points = coords_npy.shape[0];
	    unsigned int dim = coords_npy.shape[1];
	    PointList coords(num_points);
	    for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
	    coords_npy.destruct();

        cnpy::NpyArray ma_coords_in_npy = cnpy::npy_load( input_path_ma_coords_in.c_str() );
        float* ma_coords_in_carray = reinterpret_cast<float*>(ma_coords_in_npy.data);
        PointList ma_coords_in(ma_coords_in_npy.shape[0]);
        for ( int i=0; i<num_points; i++) ma_coords_in[i] = Point(&ma_coords_in_carray[i*3]);
        ma_coords_in_npy.destruct();    
        cnpy::NpyArray ma_coords_out_npy = cnpy::npy_load( input_path_ma_coords_out.c_str() );
        float* ma_coords_out_carray = reinterpret_cast<float*>(ma_coords_out_npy.data);
        PointList ma_coords_out(ma_coords_out_npy.shape[0]);
        for ( int i=0; i<num_points; i++) ma_coords_out[i] = Point(&ma_coords_out_carray[i*3]);
        ma_coords_out_npy.destruct();

        cnpy::NpyArray ma_qidx_in_npy = cnpy::npy_load( input_path_ma_qidx_in.c_str() );
        madata.ma_qidx_in = reinterpret_cast<int*>(ma_qidx_in_npy.data);
        for ( int i=0; i<num_points; i++) {
            std::cout << "ma_qidx_in: " << madata.ma_qidx_in[i] << std::endl;
        }
        // ma_qidx_in_npy.destruct

        cnpy::NpyArray ma_qidx_out_npy = cnpy::npy_load( input_path_ma_qidx_out.c_str() );
        madata.ma_qidx_out = reinterpret_cast<int*>(ma_qidx_out_npy.data);
        // ma_qidx_out_npy.destruct
	    
        
        madata.m = num_points;
        madata.coords = &coords;
        // madata.normals = &normals;
        madata.ma_coords_in = &ma_coords_in;
        madata.ma_coords_out = &ma_coords_out;

        madata.mask = new bool[num_points];
        madata.lfs = new float[num_points];

	    // omp_set_num_threads(4);

	    {
            // compute lfs
            compute_lfs(madata, cos((1.0 / 180.0) * 3.1415));

            // create filter ...
	    
	        const unsigned int c_size = num_points;
            const unsigned int shape[] = {c_size};
            cnpy::npy_save(output_filtermask.c_str(), madata.mask, shape, 1, "w");
	        cnpy::npy_save(output_lfs.c_str(), madata.lfs, shape, 1, "w");
	    }

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
