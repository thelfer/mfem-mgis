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

# Installing with SPACK

A simple way to install via spack MFEM-MGIS-MFront is the following:
~~~~{.bash}
$ git clone https://github.com/spack/spack.git
    # Please install spack outside the source directory of mfem-mgis,
    # because it can lead to some caveats using CMake
$ export SPACK_ROOT=$PWD/spack
$ source ${SPACK_ROOT}/share/spack/setup-env.sh
$ spack compiler find # Detect the available compilers on the system.
    # Please select a version that provides C and C++ and fortran compilers.
    # At this stage one can remove unwanted compiler with "spack compiler remove <XXX>".
$ spack external find m4 openssl automake ncurses
    # "external find" command tells spack machinery to detect already installed 
    # libraries/program. If such libraries/program are found they are not reinstalled.
$ spack external find autoconf libtool xz gmake cmake
$ spack external find tar tcl perl curl zlib openblas
$ cd mfem-mgis
$ git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
$ spack repo add spack-repo-mfem-mgis
$ spack install -j 8 mfem-mgis^mfem~mpi+suite-sparse
$ spack load mfem-mgis^mfem~mpi+suite-sparse
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=../install
$ make -j 4 check
$ make install
~~~~

# Creating a simple example based on `mfem-mgis`

Through the `make install` command, a simple example has been created in
your installation directory in "install/share/mfem-mgis/examples" directory.

You can copy it elsewhere together with the `env.sh` file. The example
can be compiled either using the build systems `cmake` or`make`.

## Building the example using the `cmake` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ mkdir build
$ cd build
$ cmake ..
$ make; make check
~~~~

The example may also be run as follows:

~~~~{.bash}
$ ./UniaxialTensileTest 
~~~~

You can then modify the source file and design your
own case of study.

## Building the example using the `make` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ make
$ ./UniaxialTensileTest 
~~~~

### Building in debug mode

The example and the `MFront` behaviour can be compiled in `debug` mode
by changing the call to make as follows:

~~~~{.bash}
$ make clean; make DEBUG=1
~~~~
