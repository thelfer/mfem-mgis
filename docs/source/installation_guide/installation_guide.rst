==================
Installation guide
==================

.. contents::
    :depth: 2
    :local:

Quick start
-----------

``mfem-mgis`` is packaged in `Spack <https://spack.io/>`_, which also
installs all its dependencies. These commands set up Spack:

.. code:: sh

   git clone --depth=2 --branch=v1.2.2 https://github.com/spack/spack.git
   source spack/share/spack/setup-env.sh
   spack repo update builtin --branch develop

The last command selects the ``develop`` branch of the Spack packages,
because their releases do not provide ``mfem-mgis`` yet.

The ``@master`` version of ``mfem-mgis``, its development version, is
currently recommended, as it includes many fixes:

.. code:: sh

   spack install mfem-mgis@master

The last release (1.0.4) is installed by:

.. code:: sh

   spack install mfem-mgis

Spack needs a few system tools, such as ``git``, ``python3`` and C, C++
and Fortran compilers. They are listed in the `Spack documentation
<https://spack.readthedocs.io/en/latest/installing_prerequisites.html>`_.
The first installation builds about 90 packages from their sources. It
takes about 15 minutes on a computer with 8 cores.

.. note::

   Install Spack outside the sources of ``mfem-mgis``, to avoid issues
   with CMake.

Installing with Spack
---------------------

Versions
^^^^^^^^

This documentation is generated from the ``master`` branch, so it
describes the ``@master`` version. Some of the features it describes are
not available in the last release (1.0.4). They are marked as such.

Variants
^^^^^^^^

The package has three variants:

- ``mpi`` enables MPI parallelism. It is enabled by default. Only the
  ``@master`` version can be built without MPI, with ``~mpi``.
- ``mumps`` adds the ``MUMPS`` solver to ``MFEM``. It is enabled by
  default when MPI is.
- ``int64`` uses 64 bit integers in ``hypre`` and ``metis``. It is
  enabled by default.

For example, this command installs the ``@master`` version without MPI:

.. code:: sh

   spack install mfem-mgis@master~mpi

Compilers
^^^^^^^^^

``mfem-mgis`` requires GCC 11 or newer. Spack detects the compilers of
the system the first time it runs. A compiler made available later, for
example by ``module load``, must be added with:

.. code:: sh

   spack compiler find

Using the installed package
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: sh

   spack load mfem-mgis

This command makes ``mfront`` and ``mpirun`` available. It also sets the
``MFEMMGIS_DIR`` variable, which CMake projects use to find
``mfem-mgis``. When several versions are installed, the version must be
given, for example ``spack load mfem-mgis@master``.

Building from the sources
-------------------------

To work on the sources of ``mfem-mgis``, Spack only installs its
dependencies, in an environment:

.. code:: sh

   git clone https://github.com/thelfer/mfem-mgis.git
   spack env create mfem-mgis-deps
   spack -e mfem-mgis-deps add mfem-mgis@master
   spack -e mfem-mgis-deps install --only dependencies
   spack env activate mfem-mgis-deps

Once activated, the environment gives CMake access to the dependencies.
CMake builds the sources, and the ``check`` target runs the tests:

.. code:: sh

   cmake -S mfem-mgis -B build -DCMAKE_INSTALL_PREFIX=install
   cmake --build build -j 8
   cmake --build build --target check
   cmake --build build --target install

Creating a simple example based on ``mfem-mgis``
------------------------------------------------

The installation of ``mfem-mgis`` includes a simple example, ``ex1``.
It requires the ``@master`` version. It can be copied to another location
and built with CMake:

.. code:: sh

   spack load mfem-mgis
   cp -r `spack location -i mfem-mgis`/share/mfem-mgis/examples/ex1 .
   cd ex1
   cmake -B build
   cmake --build build
   ctest --test-dir build

The example also provides a ``Makefile``, which needs the variables set
by ``share/mfem-mgis/examples/env.sh``.

Its sources can be modified to develop your own study cases.


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
   export MY_DEST=/ccc/scratch/cont002/den/${MY_LOG}/mini-test

Download Spack, mfem-mgis, and mfem-mgis-examples (not required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before proceeding, make sure to source Spack and clear your local ``~/.spack`` repository.
Download Spack and the git directories of ``mfem-mgis`` and ``mfem-mgis-examples``:

.. code-block:: bash

   cd $MY_DIR
   git clone --depth=2 --branch=v1.2.2 https://github.com/spack/spack.git
   rm -r ~/.spack
   export SPACK_ROOT=$PWD/spack
   source ${SPACK_ROOT}/share/spack/setup-env.sh
   cd ${SPACK_ROOT}
   git clone --branch=develop https://github.com/spack/spack-packages.git
   spack repo set --destination "${PWD}/spack-packages" builtin
   cd ..
   git clone https://github.com/thelfer/mfem-mgis.git
   git clone https://github.com/latug0/mfem-mgis-examples.git

Create a Spack Mirror on Your Machine (Local)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``develop`` branch of the Spack packages provides ``mfem-mgis``. The
examples require the ``@master`` version of ``mfem-mgis``. Create a
``spack`` mirror and a bootstrap directory.

.. code-block:: bash

   spack bootstrap mirror --binary-packages my_bootstrap
   spack mirror create -d re2c_mirror re2c@3.0
   cp -r re2c_mirror/_source-cache/archive/b3/ my_bootstrap/bootstrap_cache/_source-cache/archive
   spack mirror create -d mirror-mfem-mgis -D mfem-mgis@master

It’s possible that you will need some packages in your mirror, you can
specify them with the following command:

.. code-block:: bash

   spack mirror create -d mirror-mfem-mgis -D mfem-mgis@master zlib ca-certificates-mozilla zlib-ng util-macros pkgconf findutils libpciaccess libedit libxcrypt bison libevent numactl

**Copy Data to Topaze**

You’ll need to copy the following files to Topaze: 

- spack
- mfem-mgis
- mfem-mgis-examples
- mirror-mfem-mgis
- my_bootstrap

Create an archive for these files:

.. code-block:: bash

   cd $MY_DIR
   tar cvf archive.tar.gz mfem-mgis/ mfem-mgis-examples/ mirror-mfem-mgis/ spack/ my_bootstrap/
   scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$MY_DEST/

**Load Topaze modules**

Log on ``Topaze``:

.. code-block:: bash

   ssh -Y $MY_LOG@topaze.ccc.cea.fr

Load the required modules on Topaze:

.. code-block:: bash

   module load gnu/13.2.0
   module load mpi/openmpi/4.0.5
   module load cmake/3.29.6

Install mfem-mgis on Topaze
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that the installation is performed in your scratch directory, and
files are automatically removed after 3 months.

**Setup spack**

If you had a previous installation of ``spack``, please clean the environment to avoid any conflicts 

.. code-block:: bash

   rm -rf ~/.spack

Now extract the files and set-up the bootstrapping

.. code-block:: bash

   cd $MY_DEST
   tar xvf archive.tar.gz
   source $PWD/spack/share/spack/setup-env.sh
   spack repo set --destination "${PWD}/spack/spack-packages" builtin
   spack bootstrap reset -y
   spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
   spack bootstrap add --scope=site --trust local-sources $PWD/my_bootstrap/metadata/sources/
   spack buildcache update-index $PWD/my_bootstrap/bootstrap_cache
   spack bootstrap disable --scope=site github-actions-v2
   spack bootstrap disable --scope=site github-actions-v0.6
   spack bootstrap disable --scope=site spack-install
   spack bootstrap root $PWD/spack/bootstrap

Now you can look for the compilers

.. code-block:: bash

   spack compiler find

and remove the unnecessary ones ``spack`` might have found by editing the configuration file, e.g., 

.. code-block:: bash

   vim $HOME/.spack/packages.yaml

Now everything is set to bootstrap

.. code-block:: bash

   spack bootstrap now
   spack bootstrap status

If everything goes well, you will obtain something like

.. code-block:: bash

   [PASS] Core Functionalities

   [PASS] Binary packages

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

   spack mirror add MMM $PWD/mirror-mfem-mgis/

**Run installation**

.. code-block:: bash

   module load gnu/13.2.0 mpi/openmpi/4.0.5 hwloc cmake/3.29.6
   spack compiler find
   spack external find hwloc
   spack external find cmake
   spack external find openssh
   spack external find openmpi
   
   spack install mfem-mgis@master%gcc@13.2.0

Install MFEM-MGIS-examples on Topaze
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow these steps to install mfem-mgis-examples on Topaze:

.. code-block:: bash

   cd mfem-mgis-examples
   mkdir build && cd build
   spack load mfem-mgis
   cmake ..
   make -j 10
   ctest

**How to run an example (ex7)**

There are two ways to run an example, such as ex7, from its build
directory ``build/ex7``:

Using ccc_mprun
^^^^^^^^^^^^^^^

To run an example using ccc_mprun with 1024 processes and 1 core per process (-m access to other filesystem, -T time), execute the following command:

.. code-block:: bash

   ccc_mprun -n 1024 -c 1 -m work,store,scratch -T 84000 -pmilan ./mox2 -m mesh/inclusion.msh -o 1 -r 2 --post-processing 0

Using ccc_msub
^^^^^^^^^^^^^^


Here's an example of a `run.batch` job submission file to run an RVE simulation on 4096 MPI processes on the partition named milan for 86400 seconds.

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
  ccc_mprun ./mox2 -m mesh/inclusion.msh -o 1 -r 2 --post-processing 0

Then, to submit the job:

.. code-block:: bash
  
  ccc_msub run.batch


Troubleshooting
^^^^^^^^^^^^^^^

If you encounter ``Spack`` errors due to missing packages, consider the following two possibilities:

- Check if the package is already installed on Topaze by running:

.. code-block:: bash

  spack external find your-package

If the package is found, you can use it directly.

- If the package is not installed on ``Topaze``, you can add its sources to your mirror directory. If you are using an SSHFS mount, you can complete your mirror by executing the following command on your host machine:

.. code-block:: bash

  spack mirror create -d your-mirror/ -D your-package

For more questions about ``spack``, see the ``spack`` documentation.
