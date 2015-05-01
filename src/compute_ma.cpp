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
#include "kdtree2.hpp"

// typedefs
typedef float Scalar; // Scalar type for 3D points
// typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
// typedef Geometry::Vector<Scalar,3> Vector; // Type for 3D vectors
// typedef Geometry::PointTwoNTree<Point> SpatialIndex;
// typedef Geometry::ArrayKdTree<Point> SpatialIndex;
// typedef Geometry::PointKdTree<Scalar,3,Point> SpatialIndex;

// globals
const Scalar initial_radius = 1;
const Scalar delta_convergance = 1E-5;
const uint iteration_limit = 30;

inline Scalar compute_radius(Point &p, Vector &n, Point &q)
{
    // this is basic goniometry
    double d = Geometry::mag(p-q);
    Scalar cos_theta = ( n * (p-q) ) / d;
    return d/(2*cos_theta);
}

int nnn_counter =0;
double nnn_total_time =0;
inline Point sb_point(Point &p, Vector &n, SpatialIndex* kd_tree, std::vector<Point> &ma_coords)
{
    uint j=0;
    Point q(0,0,1), c;
    Scalar r, r_previous = initial_radius;
    Geometry::ClosePointSet<Point> close_points(2);

    while (1) 
    {
        // std::cout << "\nloop iteration: " << j << ", p = (" << p[0] << "," << p[1] << "," << p[2] << ") \n";
        
        // compute ball center c
        c = p - n * r_previous;

        // std::cout << "c = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";

        // find closest point to c
        Misc::Timer t1;
        kd_tree->findClosestPoints(c, close_points);
        nnn_counter++;
        q = close_points.getPoint(0);
        t1.elapse();
        nnn_total_time += t1.getTime()*1000.0;
        // std::cout<<"NN time: "<<t1.getTime()*1000.0<<" ms"<<std::endl;

        // std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";

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
                q = close_points.getPoint(1);
            }
        }
        // close_points.clear();

        // compute radius
        r = compute_radius(p,n,q);

        // std::cout << "r = " << r << "\n";

        // if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
        if (r < 0)
            r = initial_radius;
        // if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
        else if (r > initial_radius)
        {
            r = initial_radius;
            break;
        }

        // stop iteration if r has converged
        if (Math::abs(r_previous-r) < delta_convergance)
            break;

        // stop iteration if this looks like an infinite loop:
        if (j > iteration_limit)
            break;

        r_previous = r;
        j++;
        ma_coords.push_back(c);
    }

    // std::cout << j << ": (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
    return c;
}

std::vector<Point> sb_points(Point* point_array, Vector* normal_array, uint n, SpatialIndex* kd_tree)
{
    std::vector<Point> ma_coords;
    for (uint i=0; i<n; i++)
    {
        sb_point(point_array[i], normal_array[i], kd_tree, ma_coords);
        // ma_coords[i] = sb_point(point_array[i], normal_array[i], kd_tree);
    }
    std::cout << ": (#" << nnn_counter << ", " << nnn_total_time << "ms)\n";
    return ma_coords;
}

int main()
{
    cnpy::NpyArray arr = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/coords.npy");
    float* loaded_data = reinterpret_cast<float*>(arr.data);
    
    std::cout << arr.word_size << "; (" << arr.shape[0] << "," << arr.shape[1] << ")\n";

    uint num_points = arr.shape[0];
    uint dim = arr.shape[1];
    
    // Point* coords = new Point[num_points];
    // for ( int i=0; i<num_points; i++) coords[i] = Point(&loaded_data[i*3]);
    array2dfloat coords; 

    coords.resize(boost::extents[N][dim]);
    for (int i=0; i<num_points; i++) {
        for (int j=0; j<dim; j++) 
            coords[i][j] = loaded_data[i*dim+j];
    }
    // arr.destruct();
    // delete[] loaded_data;

    cnpy::NpyArray arr2 = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/normals.npy");
    float* loaded_data2 = reinterpret_cast<float*>(arr2.data);
    
    std::cout << arr2.word_size << "; (" << arr2.shape[0] << "," << arr2.shape[1] << ")\n";

    Vector* normals = new Vector[arr2.shape[0]];
    for ( int i=0; i<num_points; i++) normals[i] = Vector(&loaded_data2[i*3]);
    // arr2.destruct();
    // delete[] loaded_data2;
    

    // SpatialIndex kd_tree(num_points, coords);
    kdtree2::KDTree* tree;
    kd_tree = new kdtree2::KDTree(realdata,true);


    // Vector n(0,1,0);
    // Point p(20,100,1);
    // Point* ma_coords = new Point[num_points];
    // Point c(20,100,1);
    // Point c = sb_point(p, n, &kd_tree);
    Misc::Timer t1;
    std::vector<Point> ma_coords = sb_points(coords, normals, num_points, &kd_tree);
    t1.elapse();
    std::cout<<"NN time: "<<t1.getTime()*1000.0<<" ms"<<std::endl;
    // delete[] coords;

    Scalar* ma_coords_flat = new Scalar[num_points*3];
    for (int i=0; i<ma_coords.size(); i++)
        for (int j=0; j<3; j++)
            ma_coords_flat[i*3+j] = ma_coords[i][j];

    const unsigned int c_size = ma_coords.size();
    const unsigned int shape[] = {c_size,3};
    cnpy::npy_save("ma_coords.npy", ma_coords_flat, shape, 2, "w");
}
