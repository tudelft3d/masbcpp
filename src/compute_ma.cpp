// #include<cstdlib>
#include <iostream>
// #include <vector>
// #include <new>
// #include <map>
// #include <string>

// cnpy
#include <cnpy.h>

// Vrui
#include <Geometry/ArrayKdTree.h>
#include <Geometry/PointKdTree.h>
#include <Geometry/PointTwoNTree.h>

#include <Geometry/ClosePointSet.h>
#include <Geometry/ComponentArray.h>

#include <Math/Math.h>
#include <Misc/Timer.h>


// typedefs
typedef float Scalar; // Scalar type for 3D points
typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
typedef Geometry::Vector<Scalar,3> Vector; // Type for 3D vectors
// typedef Geometry::PointTwoNTree<Point> SpatialIndex;
typedef Geometry::ArrayKdTree<Point> SpatialIndex;
// typedef Geometry::PointKdTree<Scalar,3,Point> SpatialIndex;

// globals
const Scalar initial_radius = 1000;
const Scalar delta_convergance = 1E-5;
const uint iteration_limit = 30;

inline Scalar compute_radius(Point &p, Vector &n, Point &q)
{
    // this is basic goniometry
    double d = Geometry::mag(p-q);
    Scalar cos_theta = ( n * (p-q) ) / d;
    return d/(2*cos_theta);
}

inline Point sb_point(Point &p, Vector &n, SpatialIndex* kd_tree)
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
        // Misc::Timer t1;
        kd_tree->findClosestPoints(c, close_points);
        q = close_points.getPoint(0);
        // t1.elapse();
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
        close_points.clear();

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
    }

    return c;
}

void sb_points(Point* point_array, Vector* normal_array, uint n, SpatialIndex* kd_tree)
{
    for (uint i=0; i<n; i++)
    {
        Point c = sb_point(point_array[i], normal_array[i], kd_tree);
        // std::cout << i << ": (" << c[0] << "," << c[1] << "," << c[2] << ")\n";
    }
}

int main()
{

    cnpy::NpyArray arr = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/coords.npy");
    float* loaded_data = reinterpret_cast<float*>(arr.data);
    
    std::cout << arr.word_size << "; (" << arr.shape[0] << "," << arr.shape[1] << ")\n";

    uint num_points = 20000; //arr.shape[0];
    // std::vector<Point> pts(arr.shape[0]);
    Point* coords = new Point[num_points];
    // std::vector<Point>* coords_ptr = &coords;

    for ( int i=0; i<num_points; i++)
        // for ( int j=0; j<arr.shape[1]; j++)
        {
            coords[i] = Point(&loaded_data[i*3]);
            // std::cout << " (" << loaded_data[i*3+0] << "," << loaded_data[i*3+1] << "," << loaded_data[i*3+2] << ")\n";
            // std::cout << " (" << coords[i][0] << "," << coords[i][1] << "," << coords[i][2] << ")\n";
        }
    // arr.destruct();
    // delete[] loaded_data;

    cnpy::NpyArray arr2 = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/normals.npy");
    float* loaded_data2 = reinterpret_cast<float*>(arr2.data);
    
    std::cout << arr2.word_size << "; (" << arr2.shape[0] << "," << arr2.shape[1] << ")\n";

    // std::vector<Point> coords(arr2.shape[0]);
    Vector* normals = new Vector[arr2.shape[0]];
    // std::vector<Vector>* coords_ptr = &coords;

    for ( int i=0; i<num_points; i++)
        // for ( int j=0; j<arr2.shape[1]; j++)
        {
            normals[i] = Vector(&loaded_data2[i*3]);
            // std::cout << " (" << loaded_data2[i*3+0] << "," << loaded_data2[i*3+1] << "," << loaded_data2[i*3+2] << ")\n";
            // std::cout << " (" << coords[i][0] << "," << coords[i][1] << "," << coords[i][2] << ")\n";
        }
    // arr2.destruct();
    // delete[] loaded_data2;
    

    SpatialIndex kd_tree(num_points, coords);

    // Scalar x, y, z;
    // while (1) {
    //     std::cout << "Give me some numbers!\n";
    //     std::cin >> x >> y >> z;
    //     const Point q = Point(x,y,z);
    //     Point p = kd_tree.findClosestPoint(q);
    //     std::cout << " (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        
    // }
    

    Vector n(0,1,0);
    Point p(20,100,1);
    // Point c(20,100,1);
    // Point c = sb_point(p, n, &kd_tree);
    Misc::Timer t1;
    sb_points(coords, normals, num_points, &kd_tree);
    t1.elapse();
    std::cout<<"NN time: "<<t1.getTime()*1000.0<<" ms"<<std::endl;
    // delete[] coords;



    // cnpy::npy_save("arr_out.npy", loaded_data, &arr.shape[0], 2, "w");

    
    
}
