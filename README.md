# masbcpp
masbcpp is a C++ implementation of the shrinking ball algorithm to approximate the Medial Axis Transform (MAT) of an oriented point cloud. It is being developed in support of the [3DSM project](http://3dgeoinfo.bk.tudelft.nl/projects/3dsm/) that aims to explore possible applications of the MAT for GIS point clouds (e.g. from airborne [LiDAR](http://en.wikipedia.org/wiki/Lidar)). To deal with noisy input data a novel noise-handling mechanism is built-in. See some example output [here](https://vimeo.com/127577620).

[This video](https://vimeo.com/84859998) demonstrates how the shrinking ball algorithm works.

## Installation
For Linux/OS X:
```
$ git clone https://github.com/tudelft3d/masbcpp.git
$ cd masbcpp
$ cmake .
$ make
```
Building on windows is also possible, see the [wiki](https://github.com/tudelft3d/masbcpp/wiki/Building-on-Windows-with-CMake-and-MinGW).

###Build with OpenMP support (multithreading)
The build system should autodetect whether your compiler supports OpenMP and build masbpp with multithreading support accordingly. To enable multithreading on Mac OS X (tested with version 10.10), do
```
$ brew install clang-omp
```
prior to building masbcpp (assuming you have installed [Homebrew](http://brew.sh)).

## Usage
See
```
$ ./compute_ma --help
```
and
```
$ ./compute_normals --help
```
Currently only [NumPy](http://www.numpy.org) binary files (`.npy`) are supported as input and output. Use [lasnpyio](https://github.com/Ylannl/lasnpyio) for reading and writing of `.npy` files and conversion from the ASPRS LAS format. 

## Limitations
The current implementation is not infinitely scalable, mainly in terms of memory usage. Processing very large datasets (hundreds of millions of points or more) is therefore not really supported. 

## Acknowledgements
The shinking ball algorithm was originally introduced by [ma12]

```bib
@article{ma12,
  title={3D medial axis point approximation using nearest neighbors and the normal field},
  author={Ma, Jaehwan and Bae, Sang Won and Choi, Sunghee},
  journal={The Visual Computer},
  volume={28},
  number={1},
  pages={7--19},
  year={2012},
  publisher={Springer}
}
```

I would like to thank the authors of the following libraries for releasing their code to the general public under permissive licenses, masbcpp ships with (parts of) these excellent libraries:

* [cnpy](https://github.com/rogersce/cnpy)
* [tclap](http://tclap.sourceforge.net)
* [vrui](https://github.com/KeckCAVES/Vrui)
* [kdtree2](https://github.com/jmhodges/kdtree2)
