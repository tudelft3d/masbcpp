//
// A demonstration of using the KDTREE2 C++ routines, and timing.
// This file is in the public domain.
//

#include "kdtree2.hpp"

#include <boost/multi_array.hpp>
#include <boost/random.hpp>

#include <cnpy.h>

static boost::minstd_rand generator(42u); 
static boost::uniform_real<> uni_dist(0,1); 
boost::variate_generator<boost::minstd_rand&,boost::uniform_real<> > uni(generator,uni_dist); 

float random_variate() {
  // between [0,1)
  return(uni()); 
}

//
// define, for convenience a 2d array of floats. 
//
typedef boost::multi_array<float,2> array2dfloat;


#include <ctime>



int main() {
  // array2dfloat data(boost::extents[10][3]);  // declare a 10000 x 3 array.
  array2dfloat realdata; 

  // notice it is in C-standard layout. 
  kdtree2::KDTree* tree;
  kdtree2::KDTreeResultVector res; 
  int N, dim=3;
  int nn =2;

  /* Load or create a point set: */
  cnpy::NpyArray arr = cnpy::npy_load("/Users/ravi/git/masb/lidar/rdam_blokken_npy_lfsk10/coords.npy");
  float* loaded_data = reinterpret_cast<float*>(arr.data);    
  std::cout << arr.word_size << "; (" << arr.shape[0] << "," << arr.shape[1] << ")\n";

  cnpy::NpyArray arr2 = cnpy::npy_load("/Users/ravi/git/masbcpp/ma_coords.npy");
  float* loaded_data2 = reinterpret_cast<float*>(arr2.data);    
  std::cout << arr2.word_size << "; (" << arr2.shape[0] << "," << arr2.shape[1] << ")\n";

  N = arr.shape[0];
  realdata.resize(boost::extents[N][dim]); 
    
  for (int i=0; i<N; i++) {
    for (int j=0; j<dim; j++) 
      realdata[i][j] = loaded_data[i*dim+j];
    // std::cout <<"; (" << i << " " << N <<  " " << dim << ")\n";
  }

  int nsearch = arr2.shape[0];
  
  tree = new kdtree2::KDTree(realdata,true);
  // tree->sort_results = true;
 //  std::cout << "Tree created, now testing against brute force..."; 
 //  {
 //    std::vector<float> query(dim); 
 //    kdtree2::KDTreeResultVector result, resultbrute;
 //    int nn = 10; 

 //    for (int i=0; i<50; i++) {
 //      for (int j=0; j<dim; j++) query[j] = random_variate(); 

 //      tree->n_nearest_brute_force(query,nn,resultbrute);
 //      tree->n_nearest(query,nn,result); // search for 10 of them.

 //      for (int k=0; k<nn; k++) {
	// if ((resultbrute[k].dis != result[k].dis) ||
	//     (resultbrute[k].idx != result[k].idx)) {
	//   std::cout << "Mismatch! nn=" << k << " brute=[" << 
	//     resultbrute[k].dis << "," << resultbrute[k].idx << 
	//     "] tree=[" << result[k].dis << "," << result[k].idx << "]\n"; 
	// }
 //      }
 //    }
 //  }
  std::cout << "\nTesting complete.  Now testing timing...\n";
  tree->sort_results = false;

  {
    std::vector<float> query(dim);
    kdtree2::KDTreeResultVector result; 

    clock_t t0, t1; 

    t0 = clock();

    for (int i=0; i<nsearch;i++) {
        for (int j=0; j<dim; j++) query[j] = loaded_data2[i*dim+j]; 
        tree->n_nearest(query,nn,result);
    }

    t1 = clock(); 

    std::cout << "(" << static_cast<double> (t1-t0) / static_cast<double> (CLOCKS_PER_SEC) << ")\n";
    // return(static_cast<float> 
     // (static_cast<double> (t1-t0) / static_cast<double> (CLOCKS_PER_SEC) ));
    
  }


}

