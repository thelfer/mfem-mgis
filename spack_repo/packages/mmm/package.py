# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
# Contributions made by CEA, France

from spack import *
import os
import shutil


class Mmm(CMakePackage):
    """MMM package"""

    homepage = "https://github.com/thelfer/mfem-mgis"
    url      = "https://github.com/thelfer/mfem-mgis.git"

    version('develop', git='https://github.com/thelfer/mfem-mgis.git',
        branch='master', submodules='True')

    variant('static'      , default=True   , description='Build static library')
    variant('shared'      , default=False  , description='Build shared library')
    variant('debugmess'   , default=False  , description='Print debugging messages (add DISPMESS flag)')
    variant('exceptions'  , default=False  , description='Enable exceptions')
    variant('cuda'        , default=False  , description='Enable CUDA support')
    variant('cuda_arch'   , default='sm_60', description='CUDA architecture to compile for')
    variant('occa'        , default=False  , description='Enable OCCA backend')
    variant('raja'        , default=False  , description='Enable RAJA backend')
    variant('libceed'     , default=False  , description='Enable libCEED backend')
    variant('umpire'      , default=False  , description='Enable Umpire support')
    variant('mpi'         , default=True   , description='Enable MPI parallelism')
    variant('superlu-dist', default=False  , description='Enable MPI parallel, sparse direct solvers')
    variant('strumpack'   , default=False  , description='Enable support for STRUMPACK')
    variant('suite-sparse', default=False  , description='Enable serial, sparse direct solvers')
    variant('petsc'       , default=False  , description='Enable PETSc solvers, preconditioners, etc.')
    variant('sundials'    , default=False  , description='Enable Sundials time integrators')
    variant('pumi'        , default=False  , description='Enable functionality based on PUMI')
    variant('gslib'       , default=False  , description='Enable functionality based on GSLIB')
    variant('mpfr'        , default=False  , description='Enable precise, 1D quadrature rules')
    variant('lapack'      , default=False  , description='Use external blas/lapack routines')
    variant('debug'       , default=False  , description='Build debug instead of optimized version')
    variant('netcdf'      , default=False  , description='Enable Cubit/Genesis reader')
    variant('conduit'     , default=False  , description='Enable binary data I/O using Conduit')
    variant('zlib'        , default=True   , description='Support zip\'d streams for I/O')
    variant('gnutls'      , default=False  , description='Enable secure sockets using GnuTLS')
    variant('libunwind'   , default=False  , description='Enable backtrace on error support using Libunwind')
    variant('timer'       , default='auto' , values=('auto', 'std', 'posix', 'mac', 'mpi'), description='Timing functions to use in mfem::StopWatch')

    conflicts('~static~shared')
    conflicts('+petsc~mpi')
    conflicts('+superlu-dist~mpi')

    depends_on('mpi', when='+mpi')
    depends_on('hypre')
    depends_on('metis@5.1.0:')
    depends_on('mfem@4.7.0:')
    depends_on('mfem@4.7.0:+mpi', when='+mpi')
    depends_on('mfem@4.7.0:+suite-sparse', when='+suite-sparse')
    depends_on('mfem@4.7.0:+petsc', when='+petsc')
    depends_on('mgis@master:+c~fortran~python')
    depends_on('tfel@master:~python~python_bindings')
    depends_on('blas', when='+lapack')
    depends_on('lapack@3.0:', when='+lapack')

    depends_on('cmake', type='build')
    #    depends_on('gdb', type='build')
    depends_on('cuda', when='+cuda')
    depends_on('sundials@5.0.0:+mpi+hypre', when='+sundials')
    depends_on('sundials@5.4.0:+cuda', when='+sundials+cuda')
    depends_on('pumi@2.2.3', when='+pumi')
    depends_on('pumi', when='+pumi~shared')
    depends_on('pumi+shared', when='+pumi+shared')
    depends_on('gslib@1.0.5:+mpi', when='+gslib')
    depends_on('suite-sparse', when='+suite-sparse')
    depends_on('superlu-dist', when='+superlu-dist')
    depends_on('strumpack@3.0.0:', when='+strumpack~shared')
    depends_on('strumpack@3.0.0:+shared', when='+strumpack+shared')
    # The PETSc tests in MFEM will fail if PETSc is not configured with
    # SuiteSparse and MUMPS. On the other hand, if we require the variants
    # '+suite-sparse+mumps' of PETSc, the xsdk package concretization fails.
    depends_on('petsc@3.8:+mpi+double+hypre~hdf5', when='+petsc')
    depends_on('mpfr', when='+mpfr')
    depends_on('netcdf-c@4.1.3:', when='+netcdf')
    depends_on('unwind', when='+libunwind')
    depends_on('zlib', when='+zlib')
    depends_on('gnutls', when='+gnutls')
    depends_on('conduit@0.3.1:,master:', when='+conduit')
    depends_on('conduit+mpi', when='+conduit')


    depends_on('occa@1.1.0:', when='+occa')
    depends_on('occa+cuda', when='+occa+cuda')

    depends_on('raja@0.10.0:', when='+raja')
    depends_on('raja+cuda', when='+raja+cuda')

    depends_on('libceed@0.7:', when='+libceed')
    depends_on('libceed+cuda', when='+libceed+cuda')

    depends_on('umpire@2.0.0:', when='+umpire')
    depends_on('umpire+cuda', when='+umpire+cuda')

