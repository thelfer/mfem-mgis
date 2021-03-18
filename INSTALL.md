This project uses [`cmake`](https://cmake.org/) as build system.

# Dependencies

## Required dependencies

- [`MFEM`](https://mfem.org/)
- [`MGIS`](https://github.com/thelfer/MFrontGenericInterfaceSupport)

A simple way to install dependencies is to rely on `Spack` system. 
Spack is an open source package manager that simplifies building, installing, customizing, and sharing HPC software.
It will allow you to install recent versions of compilers (that handle c++17, for example gnu compiler suite version 8),
and to get python, cmake and other tools that are required for this project to be installed.

~~~~{.bash}
    git clone https://github.com/spack/spack.git
    export SPACK_ROOT=$PWD/spack
    source ${SPACK_ROOT}/share/spack/setup-env.sh
    spack compiler find
    spack install hypre metis mgis@master cmake
    spack load hypre metis mgis@master cmake

    git clone https://github.com/mfem/mfem.git
    # or download a tarball here : https://mfem.org/download/
    cd mfem
    mkdir build; cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=$PWD/mfem -DCMAKE_CXX_COMPILER=g++ 
    make -j 4 install
    make check
    export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
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

- Suppose that you install `mgis` using spack. For example with the command `spack install mgis@master`.
~~~~{.bash}
cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
   -DCMAKE_INSTALL_PREFIX=$PWD/../install \
   -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake

make -j 4 install
make check
~~~~

# Parallel setting

- To use parallel features of MFEM and MFEM-MGIS, you need to activate them at compile time.
  - Configuring MFEM through the following setting
~~~~{.bash}
cd mfem/build
make clean; rm CMakeCache.txt
cmake .. -DMFEM_USE_MPI=ON -DMFEM_USE_METIS_5=ON -DCMAKE_INSTALL_PREFIX=$PWD/mfem -DCMAKE_CXX_COMPILER=g++
make -j 4 install
make check
export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
~~~~
  - Configuring MFEM-MGIS with the command
~~~~{.bash}
cd mfem-mgis/build
make clean; rm CMakeCache.txt
cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
   -DCMAKE_INSTALL_PREFIX=$PWD/../install \
   -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake
make -j 4 install
make check
~~~~~

# Creating a simple app based on mfem-mgis

- Through the `make install` command, a simple application
  has been created in your installation directory.
  You can copy it elsewhere together with the `env.sh` file.
  For example :
~~~~{.bash}
export INSTALLDIR=<your_mfemmgis_install_directory>
export TGDIR=<your_target_directory>
cd ${TGDIR}
cp ${INSTALLDIR}/share/examples/env.sh .
cp -r ${INSTALLDIR}/share/examples/ex1 .
source env.sh
mkdir build
cd build
cmake ..
make
./UniaxialTensileTest 
~~~~
  You can then modify the source file and design your
  own application.

