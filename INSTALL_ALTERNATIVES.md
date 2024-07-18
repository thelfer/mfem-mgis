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
# Alternative Installation Procedure for MFEM

If you want to install and fine-tune MFEM separately, you can follow the steps outlined below. This alternative method uses Spack to manage dependencies and build MFEM with specific configurations.

## Prerequisites
- Ensure you have `git` installed on your system.

## Step 1: Clone the Spack Repository
Clone the Spack repository from GitHub.

```sh
git clone https://github.com/spack/spack.git
```

Set the `SPACK_ROOT` environment variable and source the setup script.

```sh
export SPACK_ROOT=$PWD/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

## Step 2: Detect Available Compilers
Detect the available compilers on your system.

```sh
spack compiler find
```

## Step 3: Install Required Packages
Install the required packages using Spack.

```sh
spack install hypre metis mgis@master cmake
```

## Step 4: Load the Installed Packages
Load the installed packages.

```sh
spack load hypre metis mgis@master cmake
```

## Step 5: Clone the MFEM Repository
Clone the MFEM repository from GitHub. Alternatively, you can download a tarball from [MFEM's download page](https://mfem.org/download/).

```sh
git clone https://github.com/mfem/mfem.git
```

## Step 6: Build and Install MFEM
Navigate to the `mfem` directory, create a build directory, and configure the build with CMake.

```sh
cd mfem
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PWD/mfem -DCMAKE_CXX_COMPILER=g++
```

Build and install MFEM.

```sh
make -j 4 install
```

## Step 7: Run Checks
Run the checks to verify the installation.

```sh
make check
```

## Step 8: Set MFEM_DIR Environment Variable
Set the `MFEM_DIR` environment variable to the location of the installed MFEM.

```sh
export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
```

By following these steps, you can install and fine-tune MFEM using Spack and CMake. This method allows for more control over the configuration and installation of MFEM.

# Installing MFEM-MGIS with Explicit Dependencies

This guide provides a step-by-step approach to installing MFEM-MGIS with explicit dependencies using CMake and Spack. The process involves setting relevant variables and compiling the project with specific configurations.

## Relevant Variables
- `CMAKE_BUILD_TYPE`: Type of build (e.g., Release, Debug).
- `CMAKE_INSTALL_PREFIX`: Directory where the project will be installed.
- `MFEM_DIR`: Directory where `MFEMConfig.cmake` is installed. This file is generated when MFEM is compiled and installed using CMake.
- `MFrontGenericInterface_DIR`: Directory where `MFrontGenericInterfaceConfig.cmake` is installed.

# Example of Usage

#### Step 1: Install `mgis` Using Spack
First, install `mgis` using Spack.

```sh
spack install mgis@master
```

## Step 2: Clone the MFEM Repository and Install MFEM
Follow the procedure mentioned earlier to clone, build, and install MFEM. Ensure you have the `MFEMConfig.cmake` file available.

```sh
git clone https://github.com/mfem/mfem.git
cd mfem
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PWD/mfem -DCMAKE_CXX_COMPILER=g++
make -j 4 install
make check
export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
```

## Step 3: Clone the MFEM-MGIS Repository
Clone the repository that contains MFEM-MGIS if not already done.

```sh
git clone <repository_url> mfem-mgis
cd mfem-mgis
```

## Step 4: Configure and Build MFEM-MGIS
Configure the project using CMake, specifying the required directories for MFEM and MGIS. Then, build and install the project.

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_C_COMPILER=gcc \
         -DCMAKE_CXX_COMPILER=g++ \
         -DCMAKE_INSTALL_PREFIX=$PWD/../install \
         -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake \
         -DMFEM_DIR=$MFEM_DIR
make -j 4 install
make check
```

# Explanation of the Commands

- `cmake .. -DCMAKE_BUILD_TYPE=Release ...`: Configures the build in Release mode with GCC as the compiler. The installation prefix is set to `../install`.
- `-DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake`: Specifies the directory for the `MFrontGenericInterface` configuration.
- `make -j 4 install`: Compiles and installs the project using 4 parallel jobs.
- `make check`: Runs the built-in checks to verify the installation.

By following these steps, you can explicitly manage dependencies and install MFEM-MGIS with a clear and configurable setup using Spack and CMake.
# Configuring Parallel Features for MFEM and MFEM-MGIS

To utilize the parallel features of MFEM and MFEM-MGIS, specific settings must be enabled during the compilation process. This guide outlines the steps required to configure and build these libraries with parallel support.

## Step 1: Configuring MFEM with Parallel Features

### Prerequisites
Ensure you have installed the necessary dependencies (`hypre`, `metis`, `mgis`, `cmake`) and have Spack set up.

### Instructions

1. **Navigate to the Build Directory**

   ```sh
   cd mfem/build
   ```

2. **Clean the Build Directory**

   ```sh
   make clean
   rm CMakeCache.txt
   ```

3. **Configure MFEM with MPI and METIS**

   Use CMake to configure the build with MPI and METIS support. Additionally, you can enable SuiteSparse and MUMPS for extra linear solvers.

   ```sh
   cmake .. -DMFEM_USE_MPI=ON \
            -DMFEM_USE_METIS_5=ON \
            -DCMAKE_INSTALL_PREFIX=$PWD/mfem \
            -DCMAKE_CXX_COMPILER=g++ \
            -DMFEM_USE_SUITESPARSE=ON \    # Optional
            -DMFEM_USE_MUMPS=ON            # Optional
   ```

4. **Build and Install MFEM**

   ```sh
   make -j 4 install
   ```

5. **Run Checks**

   ```sh
   make check
   ```

6. **Set the MFEM_DIR Environment Variable**

   ```sh
   export MFEM_DIR=$PWD/mfem/lib/cmake/mfem
   ```

## Step 2: Configuring MFEM-MGIS with Parallel Features

### Instructions

1. **Navigate to the Build Directory**

   ```sh
   cd mfem-mgis/build
   ```

2. **Clean the Build Directory**

   ```sh
   make clean
   rm CMakeCache.txt
   ```

3. **Configure MFEM-MGIS with CMake**

   Use CMake to configure the build, specifying the compilers and installation directories.

   ```sh
   cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_CXX_COMPILER=g++ \
            -DCMAKE_INSTALL_PREFIX=$PWD/../install \
            -DMFrontGenericInterface_DIR=$(spack location -i mgis@master)/share/mgis/cmake \
            -DMFEM_DIR=$MFEM_DIR
   ```

4. **Build and Install MFEM-MGIS**

   ```sh
   make -j 4 install
   ```

5. **Run Checks**

   ```sh
   make check
   ```

## Explanation of Additional Flags

- `-DMFEM_USE_MPI=ON`: Enables MPI support for parallel computing.
- `-DMFEM_USE_METIS_5=ON`: Enables METIS support for mesh partitioning.
- `-DMFEM_USE_SUITESPARSE=ON`: Enables SuiteSparse support for additional linear solvers (optional).
- `-DMFEM_USE_MUMPS=ON`: Enables MUMPS support for additional linear solvers (optional).

