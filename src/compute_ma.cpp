#include <iostream>

// Vrui
#include "../thirdparty/vrui/Geometry/ComponentArray.h"
#include "../thirdparty/vrui/Math/Math.h"
#include "../thirdparty/vrui/Misc/Timer.h"

// kdtree2
#include "../thirdparty/kdtree2/kdtree2.hpp"
// cnpy
#include "../thirdparty/cnpy/cnpy.h"

// typedefs
#include "types.h"

// #define VERBOSEPRINT 1;

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

double denoise_preserve = (3.1415/180) * 20;
double denoise_planar = (3.1415/180) * 32;
Point  sb_point(Point &p, Vector &n, kdtree2::KDTree* kd_tree)
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
        if (q == p)
        {
            // 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
            if (r_previous == initial_radius)
            {
                r = initial_radius;
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

    return c;
}

PointList sb_points(PointList &points, VectorList &normals, kdtree2::KDTree* kd_tree, bool inner=1)
{
    PointList ma_coords(points.size());
    Point p;
    Vector n;
    for (uint i=0; i<points.size(); i++)
    {
        p = points[i];
        if (inner)
            n = normals[i];
        else
            n = -normals[i];
        ma_coords[i] = sb_point(p, n, kd_tree);
    }
    return ma_coords;
}

int main()
{
    cnpy::NpyArray coords_npy = cnpy::npy_load("rdam_blokken_npy/coords.npy");
    float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

    uint num_points = coords_npy.shape[0];
    uint dim = coords_npy.shape[1];
    PointList coords(num_points);
    for ( int i=0; i<num_points; i++) coords[i] = Point(&coords_carray[i*3]);
    coords_npy.destruct();

    cnpy::NpyArray normals_npy = cnpy::npy_load("rdam_blokken_npy/normals.npy");
    float* normals_carray = reinterpret_cast<float*>(normals_npy.data);
    VectorList normals(normals_npy.shape[0]);
    for ( int i=0; i<num_points; i++) normals[i] = Vector(&normals_carray[i*3]);
    normals_npy.destruct();
    
    kdtree2::KDTree* kd_tree;
    kd_tree = new kdtree2::KDTree(coords,true);
    kd_tree->sort_results = true;

    {
        Scalar* ma_coords_in_carray = new Scalar[num_points*3];
        Misc::Timer t1;
        PointList ma_coords_in = sb_points(coords, normals, kd_tree, 1);
        t1.elapse();
        std::cout<<"NN time in: "<<t1.getTime()*1000.0<<" ms"<<std::endl;
    
        
        for (int i=0; i<ma_coords_in.size(); i++)
            for (int j=0; j<3; j++)
                ma_coords_in_carray[i*3+j] = ma_coords_in[i][j];
    
        const unsigned int c_size = ma_coords_in.size();
        const unsigned int shape[] = {c_size,3};
        cnpy::npy_save("rdam_blokken_npy/ma_coords_in.npy", ma_coords_in_carray, shape, 2, "w");
    }

    {
        Scalar* ma_coords_out_carray = new Scalar[num_points*3];
        Misc::Timer t2;
        PointList ma_coords_out = sb_points(coords, normals, kd_tree, 0);
        t2.elapse();
        std::cout<<"NN time out: "<<t2.getTime()*1000.0<<" ms"<<std::endl;

        
        for (int i=0; i<ma_coords_out.size(); i++)
            for (int j=0; j<3; j++)
                ma_coords_out_carray[i*3+j] = ma_coords_out[i][j];

        const unsigned int c_size = ma_coords_out.size();
        const unsigned int shape[] = {c_size,3};
        cnpy::npy_save("rdam_blokken_npy/ma_coords_out.npy", ma_coords_out_carray, shape, 2, "w");
    }

    // return 1;
}
