This project uses [`cmake`](https://cmake.org/) as build system.

# Dependencies

## Required dependencies

- [`MFEM`](https://mfem.org/)
- [`MGIS`](https://github.com/thelfer/MFrontGenericInterfaceSupport)

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
cmake .. -Denable-openmp=OFF  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PWD/../install -DMFEM_DIR=<MFEM_DIR> -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake
~~~~