---
title: Installation guide
author: Guillaume Latu, Thomas Helfer, Raphaël Prat
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

- [`MFEM`](https://mfem.org/)
- [`MGIS`](https://github.com/thelfer/MFrontGenericInterfaceSupport)

A simple way to install dependencies is to rely on [`Spack` packaging
system](https://spack.io/). `Spack` is an open source package manager
that simplifies building, installing, customizing, and sharing HPC
software. It will allow you to install recent versions of compilers
(that handle `C++17`, for example gnu compiler suite version 8), and to
get `python`, `cmake` and other tools that are required for this project
to be installed (see hereafter). Other ways to install MFEM and MGIS are
available in INSTALL_ALTERNATIVES.md.


# Installation Tutorial for `MFEM-MGIS-MFront` using Spack

This tutorial provides detailed instructions on how to install `MFEM-MGIS-MFront` using Spack. Follow the steps carefully to ensure a successful installation.

## Prerequisites
- Ensure you have `git` installed on your system.

## Step 1: Clone the Spack Repository
First, clone the Spack repository from GitHub.

```sh
git clone https://github.com/spack/spack.git
```

**Note**: Install Spack outside the source directory of `mfem-mgis` to avoid issues with CMake.

## Step 2: Set Up Spack Environment
Set up the Spack environment by configuring the `SPACK_ROOT` environment variable and sourcing the setup script.

```sh
export SPACK_ROOT=$PWD/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

## Step 3: Detect Available Compilers
Detect the available compilers on your system. Ensure to select a version that provides C, C++, and Fortran compilers.

```sh
spack compiler find
```

If necessary, remove unwanted compilers with the command:

```sh
spack compiler remove <compiler_name>
```

## Step 4: Detect External Dependencies
Use Spack to detect already installed libraries and programs to avoid reinstalling them.

```sh
spack external find m4 openssl automake ncurses
spack external find autoconf libtool xz gmake cmake
spack external find tar tcl perl curl zlib openblas
```

## Step 5: Install MFEM-MGIS-MFront
Get the `spack-repo-mfem-mgis` directory, add the Spack repository, and install the package.

```sh
git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
spack repo add spack-repo-mfem-mgis
spack install -j 8 mfem-mgis^mfem+mpi+suite-sparse
```

## Step 6: Load the Installed Package
Load the installed package.

```sh
spack load mfem-mgis^mfem~mpi+suite-sparse
```

#### Step 7: Build and Install the Project
Create a build directory, configure the project with CMake, build it, and install.

```sh
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j 4 check
make install
```

# Alternative Installation Method

If you already have `mfem`, `tfel`, and `mgis` installed via Spack, follow these steps:

## Step 1: Install Required Packages
Install the required packages using Spack.

```sh
spack install mfem+mpi+suite-sparse
spack install tfel@master:~python~python_bindings
spack install mgis@master:+c~fortran~python
```

## Step 2: Load the Installed Packages
Load the installed packages.

```sh
spack load mfem
spack load tfel
spack load mgis
spack load hypre
```

## Step 3: Set HYPRE_DIR Environment Variable
Set the `HYPRE_DIR` environment variable to the installation location of `hypre`.

```sh
export HYPRE_DIR=`spack location -i hypre`
```

## Step 4: Build and Install the Project
Create a build directory, configure the project with CMake, and build it.

```sh
mkdir build && cd build
cmake ..
make -j 4 check
```

By following these detailed instructions, you should be able to install and configure `MFEM-MGIS-MFront` using Spack successfully.


# Creating a Simple Example Based on `mfem-mgis`

Upon executing the `make install` command during the installation process, a simple example is created in your installation directory. This example can be found in the "install/share/mfem-mgis/examples" directory. You can copy this example and the associated `env.sh` file to another location. The example can be compiled using either the `cmake` or `make` build systems.

## Step 1: Locate and Copy Example Files
First, locate your installation directory and copy the example and environment setup file to a new location.

```sh
export INSTALLDIR=<your_mfemmgis_install_directory>
cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/
```

# Building the Example Using the `cmake` Build-System

## Step 1: Set Up and Compile the Example
Navigate to the example directory, source the environment setup file, create a build directory, and compile the example using `cmake`.

```sh
cd ex1
source env.sh
mkdir build
cd build
cmake ..
make
make check
```

## Step 2: Run the Example
After successfully building the example, you can run it using the following command:

```sh
./UniaxialTensileTest
```

You can then modify the source file to design your own case study.

# Building the Example Using the `make` Build-System

## Step 1: Set Up and Compile the Example

Navigate to the example directory, source the environment setup file, and compile the example using `make`.

```sh
cd ex1
source env.sh
make
```

## Step 2: Run the Example
After successfully building the example, you can run it using the following command:

```sh
./UniaxialTensileTest
```

# Building in Debug Mode

To compile the example and the `MFront` behavior in debug mode, use the following command:

```sh
make clean
make DEBUG=1
```

By following these steps, you can successfully create, build, and run a simple example based on `mfem-mgis`. Modify the source files as needed to develop and test your own study cases.
