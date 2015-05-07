// #include<cstdlib>
#include <iostream>
#include <vector>
// #include <new>
// #include <map>
// #include <string>

// cnpy
#include <cnpy.h>

// Vrui
// #include <Geometry/ArrayKdTree.h>
// #include <Geometry/PointKdTree.h>
// #include <Geometry/PointTwoNTree.h>

// #include <Geometry/ClosePointSet.h>
#include <Geometry/ComponentArray.h>

#include <Math/Math.h>
#include <Misc/Timer.h>

// kdtree2
#include <boost/multi_array.hpp>
#include "../kdtree2/kdtree2.hpp"

#include "types.h"

// #define VERBOSEPRINT 1;

// typedefs

// typedef Geometry::PointTwoNTree<Point> SpatialIndex;
// typedef Geometry::ArrayKdTree<Point> SpatialIndex;
// typedef Geometry::PointKdTree<Scalar,3,Point> SpatialIndex;

typedef boost::multi_array<float,2> array;
// typedef array::index_gen idx;
// typedef array::array_view<1>::type vec;

// globals
const Scalar initial_radius = 100;
const Scalar delta_convergance = 1E-5;
const uint iteration_limit = 30;

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
    if (result > 1) return 1;
    else if (result < -1) return -1;
    return result;
}

int nnn_counter =0;
double nnn_total_time =0;
double denoise_preserve = (3.1415/180) * 20;
double denoise_planar = (3.1415/180) * 32;
Point  sb_point(Point &p, Vector &n, kdtree2::KDTree* kd_tree, std::vector<Point> &ma_coords)
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
        Misc::Timer t1;
        kdtree2::KDTreeResultVector result;
        std::vector<Scalar> c_vec;
        for (int i=0; i<3; i++)
            c_vec.push_back(c[i]);
        // std::vector<Scalar> c_vec = {c[0], c[1], c[2]};
        kd_tree->n_nearest(c_vec,2,result);
        // kd_tree->findClosestPoints(c, close_points);
        nnn_counter++;
        // q = close_points.getPoint(0);
        for (int i=0; i<3; i++)
            q[i] = kd_tree->the_data[ result[0].idx ][i];
        t1.elapse();
        nnn_total_time += t1.getTime()*1000.0;
        // std::cout<<"NN time: "<<t1.getTime()*1000.0<<" ms"<<std::endl;

        #ifdef VERBOSEPRINT
        std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";
        #endif

        // handle case when q==p
        if (q == p)
        {
            // 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
            if (r_previous == initial_radius)
            {
                r = initial_radius;
                break;
            // 2) otherwise just pick the second closest point
            } else {
                for (int i=0; i<3; i++)
                    q[i] = kd_tree->the_data[ result[1].idx ][i];
            }
        }
        // close_points.clear();

        // compute radius
        r = compute_radius(p,n,q);

        #ifdef VERBOSEPRINT
        std::cout << "r = " << r << "\n";
        #endif

        // if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
        if (r < 0)
            r = initial_radius;
        // if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
        else if (r > initial_radius)
        {
            r = initial_radius;
            break;
        }

        // Deonoising
        // compute ball center c
        c_next = p - n * r;
        if (denoise_preserve or denoise_planar)
        {
            Scalar a = cos_angle(p-c_next, q-c_next);
            Scalar separation_angle = Math::acos(a);
            

            // std::cout << j << "\n";
            // std::cout << separation_angle << " | " << denoise_preserve << " | " << denoise_planar << ".\n";

            if ( separation_angle < denoise_preserve and j>0 and r > Geometry::mag(q-p) )
            {
                // std::cout << "denoise.\n";
                // keep previous radius:
                r=r_previous;
                break;
            }
            // if ( separation_angle < denoise_planar and j<2 )
            if ( separation_angle < denoise_planar and j==0 )
            {
                // std::cout << "planar.\n";
                r= initial_radius;
                break;
            }
        }
        // stop iteration if r has converged
        if (Math::abs(r_previous-r) < delta_convergance)
            break;

        // stop iteration if this looks like an infinite loop:
        if (j > iteration_limit)
            break;

        r_previous = r;
        c = c_next;
        j++;
        
    }

    // std::cout << j << ": (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
    return c;
}

std::vector<Point> sb_points(array &point_array, array &normal_array, uint num_points, kdtree2::KDTree* kd_tree, bool inner=1)
{
    std::vector<Point> ma_coords;
    Point p;
    Vector n;
    for (uint i=0; i<num_points; i++)
    {
            p = Point(point_array[i][0], point_array[i][1], point_array[i][2]);
            if (inner)
                n = Vector(normal_array[i][0], normal_array[i][1], normal_array[i][2]);
            else
                n = -Vector(normal_array[i][0], normal_array[i][1], normal_array[i][2]);
        // std::cout << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        ma_coords.push_back(sb_point(p, n, kd_tree, ma_coords));
        // ma_coords[i] = sb_point(point_array[i], normal_array[i], kd_tree);
    }
    // std::cout << ": (#" << nnn_counter << ", " << nnn_total_time << "ms)\n";
    return ma_coords;
}

int main()
{
    cnpy::NpyArray coords_npy = cnpy::npy_load("rdam_blokken_npy/coords.npy");
    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

    uint num_points = coords_npy.shape[0];
    uint dim = coords_npy.shape[1];
    // Point* coords = new Point[num_points];
    // for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
    array coords; 

    coords.resize(boost::extents[num_points][dim]);
    for (int i=0; i<num_points; i++) {
        for (int j=0; j<dim; j++) 
            coords[i][j] = coords_carray[i*dim+j];
    }
    // coords_npy.destruct();
    // delete[] coords_carray;

    cnpy::NpyArray normals_npy = cnpy::npy_load("rdam_blokken_npy/normals.npy");
    float* normals_carray = reinterpret_cast<float*>(normals_npy.data);
    // Vector* normals = new Vector[normals_npy.shape[0]];
    // for ( int i=0; i<num_points; i++) normals[i] = Vector(&normals_carray[i*3]);
    array normals;
    normals.resize(boost::extents[num_points][dim]);
    for (int i=0; i<num_points; i++) {
        for (int j=0; j<dim; j++) 
            normals[i][j] = normals_carray[i*dim+j];
    }
    // normals_npy.destruct();
    // delete[] normals_carray;
    

    // SpatialIndex kd_tree(num_points, coords);
    kdtree2::KDTree* kd_tree;
    kd_tree = new kdtree2::KDTree(coords,true);
    kd_tree->sort_results = true;

    // Vector n(0,1,0);
    // Point p(20,100,1);
    // Point* ma_coords = new Point[num_points];
    // Point c(20,100,1);
    // Point c = sb_point(p, n, &kd_tree);
    Scalar* ma_coords_in_carray = new Scalar[num_points*3];
    Misc::Timer t1;
    std::vector<Point> ma_coords_in = sb_points(coords, normals, num_points, kd_tree, 1);
    t1.elapse();
    std::cout<<"NN time in: "<<t1.getTime()*1000.0<<" ms"<<std::endl;
    // delete[] coords;

    
    for (int i=0; i<ma_coords_in.size(); i++)
        for (int j=0; j<3; j++)
            ma_coords_in_carray[i*3+j] = ma_coords_in[i][j];

    const unsigned int c_size = ma_coords_in.size();
    const unsigned int shape[] = {c_size,3};
    cnpy::npy_save("rdam_blokken_npy/ma_coords_in.npy", ma_coords_in_carray, shape, 2, "w");

    Scalar* ma_coords_out_carray = new Scalar[num_points*3];
    Misc::Timer t2;
    std::vector<Point> ma_coords_out = sb_points(coords, normals, num_points, kd_tree, 0);
    t2.elapse();
    std::cout<<"NN time out: "<<t2.getTime()*1000.0<<" ms"<<std::endl;
    // delete[] coords;

    
    for (int i=0; i<ma_coords_out.size(); i++)
        for (int j=0; j<3; j++)
            ma_coords_out_carray[i*3+j] = ma_coords_out[i][j];

    cnpy::npy_save("rdam_blokken_npy/ma_coords_out.npy", ma_coords_out_carray, shape, 2, "w");

    // return 0;
}
