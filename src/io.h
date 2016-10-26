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

#ifndef MASBCPP_IO_
#define MASBCPP_IO_

#include <iostream>
#include <fstream>
#include <string>

#include "madata.h"

struct io_parameters {
   bool coords;
   bool normals;
   bool ma_coords;
   bool ma_qidx;
   bool lfs;
   bool mask;
};

void npy2madata(std::string input_dir_path, ma_data &madata, io_parameters &p);
void madata2npy(std::string npy_path, ma_data &madata, io_parameters &p);

// Just a convenience function, to call when necessary.
void convertNPYtoXYZ(std::string input_dir_path);

#endif
