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

template <typename T_read, typename T_return> 
T_return read_npyarray(std::string input_file_path) {
    std::replace(input_file_path.begin(), input_file_path.end(), '\\', '/');

    cnpy::NpyArray npy_array = cnpy::npy_load( input_file_path.c_str() );
    T_read* coords_carray = reinterpret_cast<T_read*>(npy_array.data);
    T_return result(npy_array.shape[0], npy_array.shape[1]);
    for (unsigned int i=0; i<result.rows(); i++)
        for (unsigned int j=0; i<result.cols(); j++)
            result[i,j] = coords_carray[i*result.cols()+j];
    npy_array.destruct();
    return result;
}

ma_data npy2madata(std::string input_dir_path, io_parameters p) {
    madata ma_data = {};

    if (p.coords) {
        madata.coords = read_npyarray<float, ArrayX3f>(input_dir_path+"/coords.npy");
        madata.m = madata.coords.rows()
    }
    if (p.normals) {
        madata.normals = read_npyarray<float, ArrayX3f>(input_dir_path+"/normals.npy");
    }

}