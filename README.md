# masbcpp
masbcpp is a c++ implementation of the shrinking ball algorithm to approximate the Medial Axis Transform (MAT) of an oriented point cloud. It is being developed in support of the [3DSM project](http://3dgeoinfo.bk.tudelft.nl/projects/3dsm/) that aims to explore possible applications of the MAT for GIS point clouds (e.g. from airborne [LiDAR](http://en.wikipedia.org/wiki/Lidar)). To deal with noisy input data a noise-handling mechanism is built-in.

[This video](https://vimeo.com/84859998) demonstrates how the shrinking ball algorithm works using an early version of masbpy.

## Installation
```
$ git clone https://github.com/tudelft3d/masbcpp.git
$ cd masbcpp
$ cmake .
$ make
```

###Build with OpenMP support (multithreading)
This is currently only tested on Mac OS 10.10. To enable multithreading support Mac OS, do
```
brew install clang-omp
```
prior to building masbcpp.

## Usage
See
```
$ ./compute_ma
```
Currently only numpy binary arrays are supported as input and output.


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

masbcpp ships with (parts of) the following excellent libraries:

* [cnpy](https://github.com/rogersce/cnpy)
* [tclap](http://tclap.sourceforge.net)
* [vrui](https://github.com/KeckCAVES/Vrui)
* [kdtree2](https://github.com/jmhodges/kdtree2)
