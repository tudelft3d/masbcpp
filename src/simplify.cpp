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
#include <array>
#include <random>

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

    // compute bisector and filter .. rebuild kdtree .. compute lfs .. compute grid .. thin each cell
    VectorList ma_bisec_in(madata.m);
    madata.ma_bisec_in = &ma_bisec_in;
    VectorList ma_bisec_out(madata.m);
    madata.ma_bisec_out = &ma_bisec_out;
    for( int i=0; i<madata.m; i++ ){

        Vector f1_in = (*madata.coords)[i] - (*madata.ma_coords_in)[i];
        Vector f2_in = (*madata.coords)[ madata.ma_qidx_in[i] ] - (*madata.ma_coords_in)[i];

        (*madata.ma_bisec_in)[i] = (f1_in+f2_in).normalize();

        Vector f1_out = (*madata.coords)[i] - (*madata.ma_coords_out)[i];
        Vector f2_out = (*madata.coords)[ madata.ma_qidx_out[i] ] - (*madata.ma_coords_out)[i];

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
        #pragma omp parallel for private(result)
        for( unsigned int i=0; i<madata.m; i++ ){
            kd_tree.n_nearest((*madata.ma_coords_in)[i], k, result);
            madata.mask[i] = false;
            // for( int j=1; j<k; j++ ){
                float bisec_angle = acos((*madata.ma_bisec_in)[result[1].idx] * (*madata.ma_bisec_in)[i]);
                if( bisec_angle < bisec_threshold )
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

    k=1;
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
        #pragma omp parallel for private(result)
        for( unsigned int i=0; i<madata.m; i++ ){
            kd_tree.n_nearest((*madata.coords)[i], k, result);
            madata.lfs[i] = sqrt(result[0].dis);
        }
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Computed LFS "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
    }


}

inline int flatindex(int ind[], int size[]){
    return ind[0] + size[0] * (ind[1] + ind[2]*size[1]);
}

void simplify(ma_data &madata, float cellsize, float epsilon){
    #ifndef __MINGW32__
    Misc::Timer t0;
    #endif
    
    Box::Size size = madata.bbox.getSize();
    Point origin = Point(madata.bbox.min);

    int resolution[3];

    for( int i=0; i<3; i++ )
        resolution[i] = int(size[i]/cellsize)+1;
    std::cout << "grid resolution: " << resolution[0] << " " << resolution[1] << " " << resolution[2] << std::endl;

    const int ncells = resolution[0]*resolution[1]*resolution[2];
    intList* grid[ncells];
    for (int i=0; i<ncells; i++) {
        grid[i] = NULL;
    }

    int idx[3], index;
    for( int i=0; i<madata.m; i++ ){
        for( int j=0; j<3; j++){
            idx[j] = int(((*madata.coords)[i][j]-origin[j]) / cellsize);
        }
        index = flatindex(idx, resolution);
        
        if( grid[index] == NULL ){
            intList *ilist = new intList;
            grid[index] = ilist; // this should be clear later
        }
        (*grid[index]).push_back(i);
    }

    float mean_lfs, target_n, A=cellsize*cellsize;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> randu(0, 1);
    
    // parallelize?
    for( int i=0; i<ncells; i++)
        if( grid[i] != NULL ){
            int n = grid[i]->size();
            float sum=0;
            for(auto i: *grid[i])
                sum += madata.lfs[i];
            mean_lfs = sum/n;

            target_n = A/pow(epsilon*mean_lfs,2);
            for(auto i: *grid[i])
                madata.mask[i] = randu(gen) <= target_n/n;
        }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Performed grid simplification in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif


}


int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory. Estimated MAT points are written to the files 'ma_coords_in.npy' and 'ma_coords_out.npy', for interior and exterior MAT respectively.", true, "", "output dir", cmd);

        TCLAP::ValueArg<double> epsilonArg("e","preserve","denoise preserve threshold",false,0.4,"double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c","planar","denoise planar threshold",false,1,"double", cmd);
        TCLAP::ValueArg<double> bisecArg("b","bisec","bisector threshold",false,3,"double", cmd);
        
        TCLAP::SwitchArg towdimSwitch("t","twodim","use 2D grid instead of 3D grid", cmd, false);

        cmd.parse(argc,argv);
        
        float epsilon = epsilonArg.getValue();
        float cellsize = cellsizeArg.getValue();
        float bisec_threshold = (bisecArg.getValue() / 180.0) * M_PI;
        
        bool twodim = towdimSwitch.getValue();

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_path_ma_coords_in = inputArg.getValue()+"/ma_coords_in.npy";
        std::string input_path_ma_coords_out = inputArg.getValue()+"/ma_coords_out.npy";
        std::string input_path_ma_qidx_in = inputArg.getValue()+"/ma_qidx_in.npy";
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
        madata.bbox = Box(Point(&coords_carray[0]), Point(&coords_carray[0]));
	    for ( int i=0; i<num_points; i++){ 
            coords[i] = Point(&coords_carray[i*3]);
            madata.bbox.addPoint(coords[i]);
        }
	    coords_npy.destruct();
        std::cout << "bbox: " << madata.bbox.min[0] << " " << madata.bbox.min[1] << " " << madata.bbox.min[2] << " " << madata.bbox.max[0] << " " << madata.bbox.max[1] << " " << madata.bbox.max[2] << std::endl;

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
            compute_lfs(madata, bisec_threshold);
            simplify(madata, cellsize, epsilon);

            // create filter ...
	    
	        const unsigned int c_size = num_points;
            const unsigned int shape[] = {c_size};
            cnpy::npy_save(output_filtermask.c_str(), madata.mask, shape, 1, "w");
	        cnpy::npy_save(output_lfs.c_str(), madata.lfs, shape, 1, "w");
	    }

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
