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

// OpenMP
#ifdef WITH_OPENMP
    #ifdef CLANG_OMP
        #include <libiomp/omp.h>
    #else
        #include <omp.h>
    #endif
#endif

// Vrui
#include <vrui/Geometry/ComponentArray.h>
#include <vrui/Math/Math.h>
#include <vrui/Misc/Timer.h>
// kdtree2
#include <kdtree2/kdtree2.hpp>
// cnpy
#include <cnpy/cnpy.h>
// tclap
#include <tclap/CmdLine.h>

// typedefs
#include "types.h"

// globals
Scalar initial_radius;
double denoise_preserve;
double denoise_planar;
const Scalar delta_convergance = 1E-5;
const uint iteration_limit = 30;
const Point nanPoint( std::numeric_limits<Scalar>::quiet_NaN() );

inline Scalar compute_radius(Point &p, Vector &n, Point &q)
{
    // this is basic goniometry
    double d = Geometry::mag(p-q);
    Scalar cos_theta = ( n * (p-q) ) / d;
    return d/(2*cos_theta);
}

inline Scalar cos_angle(Vector p, Vector q)
{
    // Calculate the cosine of angle between vector p and q, see http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
    Scalar result = p*q / ( Geometry::mag(p) * Geometry::mag(q) );
    if( result > 1 ) return 1;
    else if( result < -1 ) return -1;
    return result;
}

Point sb_point(Point &p, Vector &n, kdtree2::KDTree* kd_tree)
{
    uint j=0;
    Scalar r, r_previous = 0;
    Point q, c_next;
    Point c = p - n * initial_radius;

    while (1) 
    {
        #ifdef VERBOSEPRINT
        std::cout << "\nloop iteration: " << j << ", p = (" << p[0] << "," << p[1] << "," << p[2] << ", n = (" << n[0] << "," << n[1] << "," << n[2] << ") \n";

        std::cout << "c = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
        #endif

        // find closest point to c
        kdtree2::KDTreeResultVector result;
        kd_tree->n_nearest(c,2,result);
        q = kd_tree->the_data[ result[0].idx ];

        #ifdef VERBOSEPRINT
        std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";
        #endif

        // handle case when q==p
        if( q == p )
        {
            // 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
            if( r_previous == initial_radius )
            {
                r = initial_radius;
                c = nanPoint;
                break;
            // 2) otherwise just pick the second closest point
            } else {
                q = kd_tree->the_data[ result[1].idx ];
            }
        }

        // compute radius
        r = compute_radius(p,n,q);

        #ifdef VERBOSEPRINT
        std::cout << "r = " << r << "\n";
        #endif

        // if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
        if( r < 0 )
            r = initial_radius;
        // if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
        else if( r > initial_radius )
        {
            r = initial_radius;
            c = nanPoint;
            break;
        }

        // compute next ball center
        c_next = p - n * r;

        // denoising
        if( denoise_preserve or denoise_planar )
        {
            Scalar a = cos_angle(p-c_next, q-c_next);
            Scalar separation_angle = Math::acos(a);

            if( denoise_preserve and ( separation_angle < denoise_preserve and j>0 and r > Geometry::mag(q-p) ) )
            {
                // keep previous radius:
                r = r_previous;
                break;
            }
            if( denoise_planar and ( separation_angle < denoise_planar and j==0 ) )
            {
                r = initial_radius;
                c = nanPoint;
                break;
            }
        }

        // stop iteration if r has converged
        if( Math::abs(r_previous-r) < delta_convergance )
            break;

        // stop iteration if this looks like an infinite loop:
        if( j > iteration_limit )
            break;

        r_previous = r;
        c = c_next;
        j++;
    }
        
    return c;
}

PointList sb_points(PointList &points, VectorList &normals, kdtree2::KDTree* kd_tree, bool inner=1)
{
    PointList ma_coords(points.size());
    Point p;
    Vector n;

    #pragma omp parallel for private(p, n)
    for( uint i=0; i<points.size(); i++ )
    {
        p = points[i];
        if( inner )
            n = normals[i];
        else
            n = -normals[i];
        ma_coords[i] = sb_point(p, n, kd_tree);
    }
    return ma_coords;
}


int main(int argc, char **argv)
{
    // parse command line arguments
    try {
        TCLAP::CmdLine cmd("Computes a MAT point approximation", ' ', "0.1");

        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file", true, "", "input dir", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output directory", true, "", "output dir", cmd);

        TCLAP::ValueArg<double> denoise_preserveArg("d","preserve","denoise preserve threshold",false,20,"double", cmd);
        TCLAP::ValueArg<double> denoise_planarArg("p","planar","denoise planar threshold",false,32,"double", cmd);
        TCLAP::ValueArg<double> initial_radiusArg("r","radius","initial ball radius",false,200,"double", cmd);

        cmd.parse(argc,argv);
        
        initial_radius = initial_radiusArg.getValue();
        denoise_preserve = (3.1415/180) * denoise_preserveArg.getValue();
        denoise_planar = (3.1415/180) * denoise_planarArg.getValue();

        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_normals_path = inputArg.getValue()+"/normals.npy";

        // check for proper in-output arguments
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
            std::string output_path = outputArg.getValue()+"/ma_coords_in.npy";
            std::ofstream outfile(output_path.c_str());    
            if(!outfile)
                throw TCLAP::ArgParseException("invalid filepath", outputArg.getValue());
        }
        
        std::cout << "Parameters: denoise_preserve="<<denoise_preserveArg.getValue()<<", denoise_planar="<<denoise_planarArg.getValue()<<", initial_radius="<<initial_radius<<"\n";

        cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
        float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

        uint num_points = coords_npy.shape[0];
        uint dim = coords_npy.shape[1];
        PointList coords(num_points);
        for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
        coords_npy.destruct();

        cnpy::NpyArray normals_npy = cnpy::npy_load( input_normals_path.c_str() );
        float* normals_carray = reinterpret_cast<float*>(normals_npy.data);
        VectorList normals(normals_npy.shape[0]);
        for ( int i=0; i<num_points; i++) normals[i] = Vector(&normals_carray[i*3]);
        normals_npy.destruct();
        
        Misc::Timer t0;
        kdtree2::KDTree* kd_tree;
        kd_tree = new kdtree2::KDTree(coords,true);
        kd_tree->sort_results = true;
        t0.elapse();
        std::cout<<"Constructed kd-tree in "<<t0.getTime()*1000.0<<" ms"<<std::endl;

        // omp_set_num_threads(4);

        {
            Misc::Timer t1;
            PointList ma_coords_in = sb_points(coords, normals, kd_tree, 1);
            t1.elapse();
            std::cout<<"Done shrinking interior balls, took "<<t1.getTime()*1000.0<<" ms"<<std::endl;
        
            Scalar* ma_coords_in_carray = new Scalar[num_points*3];   
            for (int i=0; i<ma_coords_in.size(); i++)
                for (int j=0; j<3; j++)
                    ma_coords_in_carray[i*3+j] = ma_coords_in[i][j];
        
            const unsigned int c_size = ma_coords_in.size();
            const unsigned int shape[] = {c_size,3};
            cnpy::npy_save((outputArg.getValue()+"/ma_coords_in.npy").c_str(), ma_coords_in_carray, shape, 2, "w");
        }

        {
            Misc::Timer t2;
            PointList ma_coords_out = sb_points(coords, normals, kd_tree, 0);
            t2.elapse();
            std::cout<<"Done shrinking exterior balls, took "<<t2.getTime()*1000.0<<" ms"<<std::endl;
            
            Scalar* ma_coords_out_carray = new Scalar[num_points*3];
            for (int i=0; i<ma_coords_out.size(); i++)
                for (int j=0; j<3; j++)
                    ma_coords_out_carray[i*3+j] = ma_coords_out[i][j];

            const unsigned int c_size = ma_coords_out.size();
            const unsigned int shape[] = {c_size,3};
            cnpy::npy_save((outputArg.getValue()+"/ma_coords_out.npy").c_str(), ma_coords_out_carray, shape, 2, "w");
        }

    } catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }

    return 0;
}
