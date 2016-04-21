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

// OpenMP
#ifdef WITH_OPENMP
    #include <omp.h>
#endif

// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>
#include <vrui/Geometry/PCACalculator.h>
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

Vector estimate_normal(Point &p, kdtree2::KDTree* kd_tree, int k)
{
    kdtree2::KDTreeResultVector result;
    kd_tree->n_nearest(p,k+1,result);

    Geometry::PCACalculator<3> PCACalc;
    for(int i=0; i<k+1; i++)
        PCACalc.accumulatePoint( kd_tree->the_data[ result[i].idx ] );

    double eigen_values[3];
    PCACalc.calcCovariance();
    PCACalc.calcEigenvalues(eigen_values);
    return PCACalc.calcEigenvector(eigen_values[2]);
}

VectorList estimate_normals(PointList &points, kdtree2::KDTree* kd_tree, int k)
{
    VectorList normals(points.size());

    #pragma omp parallel for
    for( unsigned int i=0; i<points.size(); i++ )
        normals[i] = estimate_normal( points[i], kd_tree, k );

    return normals;
}


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
        
        int k = kArg.getValue();

        bool kd_tree_reorder = reorder_kdtreeSwitch.getValue();

        std::string output_path = inputArg.getValue();
        if(outputArg.isSet())
            output_path = outputArg.getValue();

        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
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
        
        std::cout << "Parameters: k="<<k<<"\n";

        cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
        float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

        unsigned int num_points = coords_npy.shape[0];
        unsigned int dim = coords_npy.shape[1];
        PointList coords(num_points);
        for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
        coords_npy.destruct();

        #ifndef __MINGW32__
        Misc::Timer t0;
        #endif
        kdtree2::KDTree* kd_tree;
        kd_tree = new kdtree2::KDTree(coords,kd_tree_reorder);
        kd_tree->sort_results = false;
        #ifndef __MINGW32__
        t0.elapse();
        std::cout<<"Constructed kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;
        #endif

        // omp_set_num_threads(1);

        {
            Scalar* normals_carray = new Scalar[num_points*3];
            VectorList normals = estimate_normals(coords, kd_tree, k);
            #ifndef __MINGW32__
            t0.elapse();
            std::cout<<"Done estimating normals, took "<<t0.getTime()*1000.0<<" ms"<<std::endl;
            #endif
            
            for (int i=0; i<normals.size(); i++)
                for (int j=0; j<3; j++)
                    normals_carray[i*3+j] = normals[i][j];
        
            const unsigned int c_size = normals.size();
            const unsigned int shape[] = {c_size,3};
            cnpy::npy_save(output_path.c_str(), normals_carray, shape, 2, "w");
        }

    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
