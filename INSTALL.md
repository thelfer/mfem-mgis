---
title: Installation guide
author: Guillaume Latu, Thomas Helfer
date: 30/03/2021
lang: en-EN
link-citations: true
colorlinks: true
numbersections: true
toc: true
geometry:
  - margin=2cm
papersize: a4
figPrefixTemplate: "$$i$$"
tblPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
---

This project uses [`cmake`](https://cmake.org/) as build system.

# Dependencies

## Required dependencies

- [`MFEM`](https://mfem.org/)
- [`MGIS`](https://github.com/thelfer/MFrontGenericInterfaceSupport)

A simple way to install dependencies is to rely on [`Spack` packaging
system](https://spack.io/). `Spack` is an open source package manager
that simplifies building, installing, customizing, and sharing HPC
software. It will allow you to install recent versions of compilers
(that handle `C++17`, for example gnu compiler suite version 8), and to
get `python`, `cmake` and other tools that are required for this project
to be installed.

~~~~{.bash}
$ git clone https://github.com/spack/spack.git
$ export SPACK_ROOT=$PWD/spack
$ source ${SPACK_ROOT}/share/spack/setup-env.sh
$ spack compiler find
$ spack install hypre metis mgis@master cmake
$ spack load hypre metis mgis@master cmake
$ 
$ git clone https://github.com/mfem/mfem.git
$ # or download a tarball here : https://mfem.org/download/
$ cd mfem
$ mkdir build; cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX=$PWD/mfem -DCMAKE_CXX_COMPILER=g++ 
$ make -j 4 install
$ make check
$ export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
~~~~

## Optional dependencies

The `TFEL` project can be used for testing purposes.

# Relevant variables

- `CMAKE_BUILD_TYPE`: type of build
- `CMAKE_INSTALL_PREFIX`: installation prefix
- `MFEM_DIR` must be set to the location where `MFEMConfig.cmake` has
  been installed. To get this file, one needs to compile MFEM with CMake
  compilation process and launch a `make install` finally.
- `MFrontGenericInterface_DIR` must be set to the location where
  `MFrontGenericInterfaceConfig.cmake` has been installed.

# Example of usage

Suppose that you install `mgis` using spack. For example with the command `spack install mgis@master`.

~~~~{.bash}
$ cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
$    -DCMAKE_INSTALL_PREFIX=$PWD/../install \
$    -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake
$ 
$ make -j 4 install
$ make check
~~~~

# Parallel setting

To use parallel features of MFEM and MFEM-MGIS, you need to activate them at compile time.

1. Configuring MFEM through the following setting

~~~~{.bash}
$ cd mfem/build
$ make clean; rm CMakeCache.txt
$ cmake .. -DMFEM_USE_MPI=ON -DMFEM_USE_METIS_5=ON -DCMAKE_INSTALL_PREFIX=$PWD/mfem \
        -DCMAKE_CXX_COMPILER=g++
$ make -j 4 install
$ make check
$ export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
~~~~

   You can also add the suite-sparse and MUMPS supports for accessing two extra linear
   solvers within mfem-mgis. During the cmake configuration you just have to provide
   the flag "-DMFEM_USE_SUITESPARSE=ON" and/or "-DMFEM_USE_MUMPS=ON". 
   
2. Configuring `mfem-mgis` with the command:

~~~~{.bash}
$ cd mfem-mgis/build
$ make clean; rm CMakeCache.txt
$ cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
$       -DCMAKE_INSTALL_PREFIX=$PWD/../install \
$       -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake
$ make -j 4 install
$ make check
~~~~~

# Creating a simple example based on `mfem-mgis`

Through the `make install` command, a simple example has been created in
your installation directory.

You can copy it elsewhere together with the `env.sh` file. The example
can be compiled either using the build systems `cmake` or`make`.

## Building the example using the `cmake` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ export TGDIR=<your_work_directory>
$ cd ${TGDIR}
$ cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ mkdir build
$ cd build
$ cmake ..
$ make
~~~~

The example may then be run as follows:

~~~~{.bash}
$ ./UniaxialTensileTest 
~~~~

You can then modify the source file and design your
own case of study.

## Building the example using the `make` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ export TGDIR=<your_work_directory>
$ cd ${TGDIR}
$ cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ make
~~~~

### Building in debug mode

The example and the `MFront` behaviour may be compiled in `debug` mode
by changing the call to make as follows:

~~~~{.bash}
$ make DEBUG=1
~~~~
