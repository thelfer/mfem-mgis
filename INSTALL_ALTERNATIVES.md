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
# Other installation procedures

## MFEM

Alternatively, if one wants to install and fine tune MFEM,
it is also possible to have a different procedure given hereafter: 
A simple way to install via spack MFEM-MGIS-MFront is the following:

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

## MFEM-MGIS with explicit dependencies

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

## Parallel setting

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
