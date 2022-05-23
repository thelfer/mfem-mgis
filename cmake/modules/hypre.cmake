# Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Defines the following variables:
#   - HYPRE_FOUND
#   - HYPRE_LIBRARIES
#   - HYPRE_INCLUDE_DIRS
#   - HYPRE_VERSION
#   - HYPRE_USING_HIP (internal)

if(MFEM_USE_MPI)
  set(HYPRE_INCLUDE_DIRS "$ENV{HYPRE_DIR}/include")
  try_run(HYPRE_VERSION_RUN_RESULT HYPRE_VERSION_COMPILE_RESULT
          ${CMAKE_CURRENT_BINARY_DIR}/cmake/modules/
          ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/get_hypre_version.cpp
          CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${HYPRE_INCLUDE_DIRS}
          RUN_OUTPUT_VARIABLE HYPRE_VERSION_OUTPUT)
  if ((HYPRE_VERSION_RUN_RESULT EQUAL 0) AND HYPRE_VERSION_OUTPUT)
    string(STRIP "${HYPRE_VERSION_OUTPUT}" HYPRE_VERSION)
    set(HYPRE_VERSION ${HYPRE_VERSION} CACHE STRING "HYPRE version." FORCE)
    message(STATUS "Found HYPRE version ${HYPRE_VERSION}")
  else()
    message(FATAL_ERROR "Unable to determine HYPRE version.")
  endif()
endif(MFEM_USE_MPI)