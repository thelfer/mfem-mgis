==================
Installation guide
==================

This project uses ``cmake`` as build system.

Dependencies
------------

- `MFEM <https://mfem.org/>`_
-  `MGIS <https://github.com/thelfer/MFrontGenericInterfaceSupport>`_

A simple way to install dependencies is to rely on ``Spack`` packaging
system. ``Spack`` is an open source package
manager that simplifies building, installing, customizing, and sharing
HPC software. It will allow you to install recent versions of compilers
(that handle ``C++17``, for example gnu compiler suite version 8), and
to get ``python``, ``cmake`` and other tools that are required for this
project to be installed (see hereafter).

Other ways to install MFEM and MGIS are available in the file
``INSTALL_ALTERNATIVES.md``.

Installation Tutorial for ``MFEM-MGIS-MFront`` using Spack
----------------------------------------------------------

This tutorial provides detailed instructions on how to install
``MFEM-MGIS-MFront`` using `Spack <https://spack.io/>`_. Follow the
steps carefully to ensure a successful installation.

Prerequisites
^^^^^^^^^^^^^

- Ensure you have ``git``, ``mpi``, ``hypre``, ``tfel``, ``mgis``, and ``mfem`` installed on your system with a recent c++ gnu compiler.
- Use ``spack`` to install missing Prerequisites.

Step 1: Clone the Spack Repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, clone the Spack repository from GitHub.

.. code:: sh

   git clone --depth=2 --branch=v0.23.0 https://github.com/spack/spack.git

.. note::

  Install Spack outside the source directory of ``mfem-mgis`` to avoid issues with CMake.

Step 2: Set Up Spack Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set up the Spack environment by configuring the ``SPACK_ROOT``
environment variable and sourcing the setup script.

.. code:: sh

   export SPACK_ROOT=$PWD/spack
   source ${SPACK_ROOT}/share/spack/setup-env.sh

Step 3: Detect Available Compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Detect the available compilers on your system. Ensure to select a
version that provides C, C++, and Fortran compilers.

.. code:: sh

   spack compiler find

If necessary, remove unwanted compilers with the command:

.. code:: sh

   spack compiler remove <compiler_name>

Step 4: Detect External Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Spack to detect already installed libraries and programs to avoid
reinstalling them.

.. code:: sh

   spack external find m4 openssl automake ncurses
   spack external find autoconf libtool xz gmake cmake
   spack external find tar tcl perl curl zlib openblas

Step 5: Install MFEM-MGIS-MFront
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Change to the ``mfem-mgis`` directory, add the Spack repository, and
install the package.

.. code:: sh

   git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
   spack repo add spack-repo-mfem-mgis
   spack install -j 8 mfem-mgis

Step 6: Load the Installed Package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Load the installed package.

.. code:: sh

   spack load mfem-mgis

Step 7: Build and Install the Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a build directory, configure the project with CMake, build it,
and install.

.. code:: sh

   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=../install
   make -j 4 check
   make install

Alternative Installation Method
-------------------------------

If you already have ``mfem``, ``tfel``, and ``mgis`` installed via
Spack, follow these steps:

Step 1: Install Required Packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the required packages using Spack.

.. code:: sh

   spack install mfem+mpi+suite-sparse
   spack install tfel@master:~python~python_bindings
   spack install mgis@master:+c~fortran~python

Step 2: Load the Installed Packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Load the installed packages.

.. code:: sh

   spack load mfem
   spack load tfel
   spack load mgis
   spack load hypre

Step 3: Set HYPRE_DIR Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set the ``HYPRE_DIR`` environment variable to the installation location
of ``hypre``.

.. code:: sh

   export HYPRE_DIR=`spack location -i hypre`

Step 4: Build and Install the Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a build directory, configure the project with CMake, and build
it.

.. code:: sh

   mkdir build && cd build
   cmake ..
   make -j 4 check

By following these detailed instructions, you should be able to install
and configure ``MFEM-MGIS-MFront`` using Spack successfully.

Creating a Simple Example Based on ``mfem-mgis``
------------------------------------------------

Upon executing the ``make install`` command during the installation
process, a simple example is created in your installation directory.
This example can be found in the “install/share/mfem-mgis/examples”
directory. You can copy this example and the associated ``env.sh`` file
to another location. The example can be compiled using either the
``cmake`` or ``make`` build systems.

Step 1: Locate and Copy Example Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, locate your installation directory and copy the example and
environment setup file to a new location.

.. code:: sh

   export INSTALLDIR=<your_mfemmgis_install_directory>
   cp -r ${INSTALLDIR}/share/mfem-mgis/examples/ex1 .
   cp ${INSTALLDIR}/share/mfem-mgis/examples/env.sh ex1/

Step 2: Set Up and Compile the Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Navigate to the example directory, source the environment setup file,
create a build directory, and compile the example using ``cmake``.

.. code:: sh

   cd ex1
   source env.sh
   mkdir build
   cd build
   cmake ..
   make
   make check

Step 3: Run the Example
^^^^^^^^^^^^^^^^^^^^^^^

After successfully building the example, you can run it using the
following command:

.. code:: sh

   ./UniaxialTensileTest

Building in Debug Mode
----------------------

To compile the example and the ``MFront`` behavior in debug mode, use
the following command:

.. code:: sh

   make clean
   make DEBUG=1

By following these steps, you can successfully create, build, and run a
simple example based on ``mfem-mgis``. Modify the source files as needed
to develop and test your own study cases.



Installation Guide on Topaze/CCRT of mfem-mgis-examples
-------------------------------------------------------

This guide provides step-by-step instructions for setting up your
environment on ``Topaze/CCRT`` and installing the necessary software. Follow
these steps to get started.

Create a new directory and useful paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   mkdir topaze-dir && cd topaze-dir
   export MY_DIR=$PWD
   export MY_LOG=YOURLOGIN
   export MY_DEST=/ccc/scratch/cont002/den/YOURLOGIN/mini-test

Download Spack, mfem-mgis, and mfem-mgis-examples (not required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to download Spack:

.. code-block:: bash

   cd $MY_DIR
   git clone --depth=2 --branch=v0.23.0 https://github.com/spack/spack.git
   export SPACK_ROOT=$PWD/spack
   git clone https://github.com/thelfer/mfem-mgis.git
   git clone https://github.com/latug0/mfem-mgis-examples.git

Before proceeding, make sure to source Spack and clear your local ``~/.spack`` repository (warning).


.. code-block:: bash

   rm -r ~/.spack
   source ${SPACK_ROOT}/share/spack/setup-env.sh

Create a Spack Mirror on Your Machine (Local)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Firstly, you need to get the mfem-mgis spack repository.

.. code-block:: bash

   git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
   spack repo add $PWD/spack-repo-mfem-mgis

Now, you will create a ``spack`` mirror and a boostrap directory. 

.. code-block:: bash

   spack bootstrap mirror --binary-packages my_bootstrap
   spack mirror create -d mirror-mfem-mgis -D mfem-mgis

It’s possible that you will need some packages in your mirror, you can
specify them with the following command:

.. code-block:: bash

   spack mirror create -d mirror-mfem-mgis -D mfem-mgis zlib ca-certificates-mozilla zlib-ng util-macros pkgconf findutils libpciaccess libedit libxcrypt bison libevent numactl

**Copy Data to Topaze**

You’ll need to copy the following files to Topaze: 

- spack spack 
- mfem-mgis 
- mfem-mgis-example

Create an archive for these files:

.. code-block:: bash

   cd $MY_DIR
   tar cvf archive.tar.gz mfem-mgis/ mfem-mgis-examples/ mirror-mfem-mgis/ spack/ my_bootstrap/ spack-repo-mfem-mgis/
   scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$MY_DEST/

**Load Topaze modules**

Log on ``Topaze``:

.. code-block:: bash

   ssh -Y $MY_LOG@topaze.ccc.cea.fr

Load the required modules on Topaze:

.. code-block:: bash

   module load gnu/11.1.0
   module load mpi

Install mfem-mgis on Topaze
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that the installation is performed in your scratch directory, and
files are automatically removed after 3 months.

**Setup spack**

.. code-block:: bash

   cd $MY_DEST
   tar xvf archive.tar.gz
   source $PWD/spack/share/spack/setup-env.sh
   spack bootstrap reset -y
   spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
   spack bootstrap disable --scope=site github-actions-v0.5
   spack bootstrap disable --scope=site github-actions-v0.6
   spack bootstrap disable --scope=site spack-install
   spack bootstrap root $PWD/spack/bootstrap
   spack repo add spack-repo-mfem-mgis/
   spack bootstrap now
   spack bootstrap status

**Export SPACK Variables**

To use ``MFront``, you need to export some ``SPACK`` variables. Please execute
the following commands:

.. code-block:: bash

   export CC='gcc'
   export CXX='g++'
   export FC='mpifort'
   export OMPI_CC='gcc'
   export OMPI_CXX='g++'
   export OMPI_FC='gfortran'

**Install MFEM-MGIS**

.. code-block:: bash

   spack repo add $PWD/spack-repo-mfem-mgis
   spack mirror add MMM $PWD/mirror-mfem-mgis/

**Run installation**

.. code-block:: bash

   module load gnu/12.3.0 mpi hwloc cmake
   spack compiler find
   spack external find hwloc
   spack external find cmake
   spack external find openssh
   spack external find openmpi
   spack install mfem-mgis%gcc@12.3.0

Install MFEM-MGIS-example on Topaze
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow these steps to install mfem-mgis-example on Topaze:

.. code-block:: bash

   cd mfem-mgis-example
   mkdir build && cd build
   spack load mfem-mgis
   export MFEMMGIS_DIR=`spack location -i mfem-mgis`/share/mfem-mgis/cmake/
   cmake ..
   make -j 10
   ctest

**How to run an example (ex8)**

There are two ways to run an example, such as ex8:

Using ccc_mprun
^^^^^^^^^^^^^^^

To run an example using ccc_mprun with 1024 processes and 1 core per process (-m access to other filesystem, -T time), execute the following command:

.. code-block:: bash

   ccc_mprun -n 1024 -c 1 -m work,store,scratch -T 84000 -pmilan ./uniaxial-elastic

Using ccc_msub
^^^^^^^^^^^^^^


Here's an example of a `run.batch` job submission file to run a VER simulation on 4096 MPI processes on the partition named milan for 84000 seconds. 

.. code-block:: bash

  #!/bin/bash
  #MSUB -r ver
  #MSUB -n 4096
  #MSUB -c 1
  #MSUB -T 86400
  #MSUB -o ver_4096_%I.o
  #MSUB -e ver_4096_%I.e
  #MSUB -q milan
  #MSUB -m scratch,work

  module load gnu/13.2.0 mpi/openmpi/4.0.5 cmake/3.29.6
  export OMP_NUM_THREADS=1
  set -x
  ccc_mprun ./ex7 -m ../par-mesh/mesh-4096. -o 1 -r 2 --post-processing 0

Then, to submit the job:

.. code-block:: bash
  
  ccc_msub run.batch


Troubleshooting
^^^^^^^^^^^^^^^

If you encounter ``Spack`` errors due to missing packages, consider the following possibilities:

Two possibilities:

- Check if the package is already installed on Topaze by running:**

.. code-block:: bash

  spack external find your-package

If the package is found, you can use it directly.

- If the package is not installed on ``Topaze``, you can add its sources to your mirror directory. If you are using an SSHFS mount, you can complete your mirror by executing the following command on your host machine:

.. code-block:: bash

  spack mirror create -d your-mirror/ -D your-package

For more questions about ``spack``, see the ``spack`` documentation.
