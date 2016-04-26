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


void compute_lfs(ma_data &madata, float bisec_threshold, bool only_inner=true)
{
    #ifndef __MINGW32__
    Misc::Timer t0;
    #endif

    int N = 2*madata.m;
    if( only_inner ){
        N = madata.m;
        (*madata.ma_coords).resize(N); // HACK this will destroy permanently the exterior ma_coords!
    }
    // compute bisector and filter .. rebuild kdtree .. compute lfs .. compute grid .. thin each cell

    VectorList ma_bisec(N);
    madata.ma_bisec = &ma_bisec;
    for( int i=0; i<N; i++ ){

        if( madata.ma_qidx[i] != -1 ){
            Vector f1_in = (*madata.coords)[ i%madata.m ] - (*madata.ma_coords)[i];
            Vector f2_in = (*madata.coords)[ madata.ma_qidx[i] ] - (*madata.ma_coords)[i];
    
            (*madata.ma_bisec)[i] = (f1_in+f2_in).normalize();
        }
        // madata.ma_theta_in = np.arccos(np.sum(f1_in*f2_in,axis=1))
    }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Computed bisectors in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif

    int k = 2, count=0;
    {
        kdtree2::KDTree kd_tree(*madata.ma_coords,true);
        kd_tree.sort_results = true;
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Constructed kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
        
        kdtree2::KDTreeResultVector result;
        #pragma omp parallel for private(result)
        for( unsigned int i=0; i<N; i++ ){
            madata.mask[i] = false;
            if( madata.ma_qidx[i] != -1 ){
                kd_tree.n_nearest((*madata.ma_coords)[i], k, result);
                
                float bisec_angle = acos((*madata.ma_bisec)[result[1].idx] * (*madata.ma_bisec)[i]);
                if( bisec_angle < bisec_threshold )
                    madata.mask[i] = true;
            }
        }
        for( unsigned int i=0; i<N; i++ )
            if (madata.mask[i])
                count++;

        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Cleaned MA points in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif
    }
    // mask and copy pointlist ma_coords
    PointList ma_coords_masked(count);
    int j=0;
    for( int i=0; i<N; i++ ){
        if( madata.mask[i] )
            ma_coords_masked[j++] = (*madata.ma_coords)[i];
    }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Copied cleaned MA points in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif

    k=1;
    {
        // rebuild kd-tree
        kdtree2::KDTree kd_tree(ma_coords_masked,true);
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

inline int flatindex(int ind[], int size[], int dimension){
    if( dimension==2 )
        return ind[0] + size[0] * ind[1];
    return ind[0] + size[0] * (ind[1] + ind[2]*size[1]);
}

void simplify(ma_data &madata, float cellsize, float epsilon, int dimension=3, float elevation_threshold=0){
    #ifndef __MINGW32__
    Misc::Timer t0;
    #endif
    
    Box::Size size = madata.bbox.getSize();
    Point origin = Point(madata.bbox.min);

    int resolution[dimension];

    for( int i=0; i<dimension; i++ )
        resolution[i] = int(size[i]/cellsize)+1;
    std::cout << "grid resolution: " << resolution[0] << " " << resolution[1] << " " << resolution[2] << std::endl;

    int ncells = 1;
    for( int i=0; i<dimension; i++ )
        ncells *= resolution[i];

    intList** grid = new intList*[ncells];
    for (int i=0; i<ncells; i++) {
        grid[i] = NULL;
    }

    int idx[dimension], index;
    for( int i=0; i<madata.m; i++ ){
        for( int j=0; j<dimension; j++){
            idx[j] = int(((*madata.coords)[i][j]-origin[j]) / cellsize);
        }
        index = flatindex(idx, resolution, dimension);
        
        if( grid[index] == NULL ){
            intList *ilist = new intList;
            // std::unique_ptr<intList> ilist(new intList);
            grid[index] = ilist;
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
            float sum=0, max_z, min_z;
            max_z = min_z = (*madata.coords)[(*grid[i])[0]][2];

            for(auto j: *grid[i]){
                sum += madata.lfs[j];
                float z = (*madata.coords)[j][2];
                if (z>max_z) max_z=z;
                if (z<min_z) min_z=z;
            }

            mean_lfs = sum/n;

            if( elevation_threshold != 0 and (max_z-min_z) > elevation_threshold )
                mean_lfs /= 5;
                

            target_n = A/pow(epsilon*mean_lfs,2);
            for(auto i: *grid[i])
                madata.mask[i] = randu(gen) <= target_n/n;
        }
    #ifndef __MINGW32__
    t0.elapse();
    std::cout<<"Performed grid simplification in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
    #endif


    // clear some memory in non-smart ptr way
    for (int i=0; i<ncells; i++) {
        delete grid[i];
    }
    delete[] grid;

}


int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to input directory with inside it a 'coords.npy' and 'ma_*.npy' files. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory", false, "", "output dir", cmd);

        TCLAP::ValueArg<double> epsilonArg("e","epsilon","espilon threshold that controls the degree of simplification",false,0.1,"double", cmd);
        TCLAP::ValueArg<double> cellsizeArg("c","cellsize","cellsize used during grid-based lfs simplification",false,1,"double", cmd);
        TCLAP::ValueArg<double> bisecArg("b","bisec","bisector threshold used to clean the MAT points",false,1,"double", cmd);
        
        TCLAP::ValueArg<double> fake3dArg("f","fake3d","use 2D grid instead of 3D grid, intended for 2.5D datasets, provide the elevation_threshold",false,0,"double", cmd);
        
        TCLAP::ValueArg<std::string> outputXYZArg("a","xyz","output filtered points to plain .xyz text file",false,"lfs_simp.xyz","string", cmd);

        cmd.parse(argc,argv);
        
        float epsilon = epsilonArg.getValue();
        float cellsize = cellsizeArg.getValue();
        float bisec_threshold = (bisecArg.getValue() / 180.0) * M_PI;
        
        float elevation_threshold = fake3dArg.getValue();
        int dimension = 3;
        bool only_inner = fake3dArg.isSet();
        if( only_inner )
            dimension = 2;

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();

        // check for proper in-output arguments and set in and output filepath strings
        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_path_ma_coords_in = inputArg.getValue()+"/ma_coords_in.npy";
        std::string input_path_ma_coords_out = inputArg.getValue()+"/ma_coords_out.npy";
        std::string input_path_ma_qidx_in = inputArg.getValue()+"/ma_qidx_in.npy";
        std::string input_path_ma_qidx_out = inputArg.getValue()+"/ma_qidx_out.npy";
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
        std::cout << "bbox: " << madata.bbox.min[0] << " " << madata.bbox.min[1] << " " << madata.bbox.min[2] << " " << madata.bbox.max[0] << " " << madata.bbox.max[1] << " " << madata.bbox.max[2] << std::endl;

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
	    
        
        madata.coords = &coords;
        // madata.normals = &normals;
        madata.ma_coords = &ma_coords;

        madata.mask = new bool[madata.m];
        madata.lfs = new float[madata.m];

	    {
            // compute lfs, simplify
            compute_lfs(madata, bisec_threshold, only_inner);
            simplify(madata, cellsize, epsilon, dimension, elevation_threshold);
	    
	        const unsigned int c_size = madata.m;
            const unsigned int shape[] = {c_size};
            cnpy::npy_save(output_filtermask.c_str(), madata.mask, shape, 1, "w");
	        cnpy::npy_save(output_lfs.c_str(), madata.lfs, shape, 1, "w");
	    }

        if( outputXYZArg.isSet() ){
            std::string outFile_bounds = outputXYZArg.getValue();
            outFile_bounds.append(".bounds");
            
            std::ofstream ofs_bounds(outFile_bounds.c_str());
            ofs_bounds << madata.bbox.min[0] << std::endl << madata.bbox.max[0] << 
             std::endl << madata.bbox.min[1] << std::endl << madata.bbox.max[1] << 
             std::endl << madata.bbox.min[2] << std::endl << madata.bbox.max[2] << std::endl;
            
            ofs_bounds.close();
            
            std::ofstream ofs(outputXYZArg.getValue());
            ofs <<std::setprecision(2)<<std::fixed;
            
            for( int i=0; i<madata.m; i++ ) {
                if( madata.mask[i] ){
                    ofs << " " << (*madata.coords)[i][0];
                    ofs << " " << (*madata.coords)[i][1];
                    ofs << " " << (*madata.coords)[i][2];
                ofs << std::endl;
                }
            }
            
            ofs.close();
        }

	} catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
