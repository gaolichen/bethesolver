# Twisted Bethe equations solver
A command line tool for finding roots of twisted Bethe equations using [Homotopy continuation](https://en.wikipedia.org/wiki/Numerical_algebraic_geometry) algorithm.

## Features
- Finding roots of Twisted Bethe equations.

## Build & Run
Recommand to run the following in linux system.

- Install [CMake](https://cmake.org/)
- Install C++ library [BOOST](https://www.boost.org/)
- Install C++ library [Eigen](https://eigen.tuxfamily.org/). Modify the CMakeList.txt in the project root directory.
```shell
SET(BOOST_ROOT ${MYDEV}/boost_1_70_0) ==> SET(BOOST_ROOT path/to/boost/library)

# If eigen directory is at the same level of the project root, then the following change is not needed.
SET(EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../eigen/) ==> SET(EIGEN_DIR path/to/eigen/library)
```
- Create a directory `build` under the project root.
- From the `build` directory, run the command line `cmake ../` then run `make`
- Run `.src/bs` for the help of the command line.

## Usage & Examples
Usage:
```shell
bs -L <spin chain length> -M <magnon number> [-b <value of beta> -r <index of root> -l <log level> -s]
```

`-s` flag indicates saving the result.

Example:
```shell
$ bs -L 8 -M 4 -b 0.3 -r 3
processing root #3:	[(0.22056,0) (0.224281,1.00225) (0.224281,-1.00225) (-0.669123,0)]
root at beta=0.3:	[(-0.509341,0) (-2.24907,2.4183) (-2.24907,-2.4183) (-2.87738,0)]
```

