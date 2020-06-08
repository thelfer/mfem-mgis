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
  been installed.
- `MFrontGenericInterface_DIR` must be set to the location where
  `MFrontGenericInterfaceConfig.cmake` has been installed.

# Example of usage

~~~~{.bash}
cmake ../mfem-mgis/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/th202608/codes/mfem-mgis/master/install -DMFEM_DIR=/home/th202608/codes/mfem/master/install/lib/cmake/mfem/ -DMFrontGenericInterface_DIR=/home/th202608/codes/mgis/master/install-python-3.5/share/mgis/cmake/
~~~~