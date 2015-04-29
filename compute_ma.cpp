// #include<cstdlib>
#include <iostream>
#include <vector>
// #include<map>
// #include<string>

// cnpy
#include <cnpy.h>

// Vrui
#include <Geometry/ArrayKdTree.h>
#include <Geometry/ClosePointSet.h>
#include <Math/Math.h>
#include <Geometry/ComponentArray.h>


// typedefs
typedef float Scalar; // Scalar type for 3D points
typedef Geometry::Point<Scalar,3> Point; // Type for 3D points
typedef Geometry::Vector<Scalar,3> Vector; // Type for 3D vectors

// globals
Geometry::ArrayKdTree<Point> kd_tree;
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

Point sb_point(Point &p, Vector &n)
{

    uint j=0;
    Point q(0,0,1), c;
    Scalar r, r_previous = initial_radius;

    Geometry::ClosePointSet<Point> close_points(2);
    // r = compute_radius(p,n,q);
    // std::cout << "radius: " << r << "\n";

    while (1) 
    {
        std::cout << "\nloop iteration: " << j << ", p = (" << p[0] << "," << p[1] << "," << p[2] << ") \n";
        if (j>0)
            r_previous = r;

        c = p - n * r_previous;

        std::cout << "c = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";

        close_points.clear();
        kd_tree.findClosestPoints(c, close_points);
        q = close_points.getPoint(0);

        std::cout << "q = (" << q[0] << "," << q[1] << "," << q[2] << ")\n";

        if (q == p)
            if (r_previous == initial_radius)
            {
                r = initial_radius;
                break;
            }
            else
            {
                q = close_points.getPoint(1);
            }

        r = compute_radius(p,n,q);

        std::cout << "r = " << r << "\n";

        if (r < 0)
            r = initial_radius;
        else if (r > initial_radius)
        {
            r = initial_radius;
            break;
        }

        if (Math::abs(r_previous-r) < delta_convergance)
            break;

        if (j > iteration_limit)
            break;

        j++;
    }

    return c;
}

int main()
{

    cnpy::NpyArray arr = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/coords.npy");
    float* loaded_data = reinterpret_cast<float*>(arr.data);
    
    std::cout << arr.word_size << "; (" << arr.shape[0] << "," << arr.shape[1] << ")\n";

    std::vector<Point> pts(arr.shape[0]);

    for ( int i=0; i<arr.shape[0]; i++)
        // for ( int j=0; j<arr.shape[1]; j++)
        {
            pts[i] = Point(&loaded_data[i*3]);
            // std::cout << " (" << loaded_data[i*3+0] << "," << loaded_data[i*3+1] << "," << loaded_data[i*3+2] << ")\n";
            // std::cout << " (" << pts[i][0] << "," << pts[i][1] << "," << pts[i][2] << ")\n";
        }
    

    kd_tree = Geometry::ArrayKdTree<Point> (arr.shape[0], &pts[0]);


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
    Point c = sb_point(p, n);
    std::cout << " (" << c[0] << "," << c[1] << "," << c[2] << ")\n";


    // cnpy::npy_save("arr_out.npy", loaded_data, &arr.shape[0], 2, "w");

    // delete[] loaded_data;
    
}
