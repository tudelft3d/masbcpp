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

#include <iostream>
#include <fstream>
#include <string>

// cnpy
#include <cnpy/cnpy.h>

// typedefs
#include "madata.h"

#include "io.h"

inline cnpy::NpyArray read_npyarray(std::string input_file_path) {
    // windows fix
    std::replace(input_file_path.begin(), input_file_path.end(), '\\', '/');
    // check if file exists
    {
        std::ifstream infile(input_file_path.c_str());    
        if(!infile){
            std::cerr << "Invalid file path " << input_file_path << std::endl;
            exit(1);
        }
    }
    // std::cout << "Reading array from " << input_file_path <<std::endl;
    cnpy::NpyArray npy_array = cnpy::npy_load( input_file_path.c_str() );
    return npy_array;
}

void npy2madata(std::string input_dir_path, ma_data &madata, io_parameters &p) {

    std::cout << "Reading coords array..." <<std::endl;
    cnpy::NpyArray npy_array = read_npyarray(input_dir_path+"/coords.npy");
    float* coords_carray = reinterpret_cast<float*>(npy_array.data);
    size_t m = npy_array.shape[0];
    for (size_t i=0; i<m; i++)
        madata.coords->push_back(Point(
            coords_carray[i*3+0],
            coords_carray[i*3+1],
            coords_carray[i*3+2]
        ));
    npy_array.destruct();

    madata.m = m;

    // if (p.normals) {
    //     std::cout << "Reading normals array..." <<std::endl;
    //     madata.normals = read_npyarray<float, ArrayX3>(input_dir_path+"/normals.npy");
    // }
}

void madata2npy(std::string npy_path, ma_data &madata, io_parameters &p) {
    if (p.normals) {
        std::cout << "Writing normals array..." <<std::endl;

        const unsigned int c_size = (unsigned int) madata.m;
        const unsigned int shape[] = { c_size, 3 };
              // Output results
        float* normals_carray = new float[c_size * 3];
        for (size_t i = 0; i < c_size; i++){
            normals_carray[i * 3 + 0] = madata.normals->points[i].normal_x;
            normals_carray[i * 3 + 1] = madata.normals->points[i].normal_y;
            normals_carray[i * 3 + 2] = madata.normals->points[i].normal_z;
        }
        cnpy::npy_save(npy_path+"/normals.npy", normals_carray, shape, 2, "w");
     }
}
