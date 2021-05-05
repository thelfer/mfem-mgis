cmake_minimum_required(VERSION 3.14)
project(mfem-mgis)

include(cmake/modules/mfem-mgis.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# (must be placed *before* any add_subdirectory, cmake bug ?)
include(CTest)
enable_testing()

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

#compiler options
include(cmake/modules/compiler.cmake)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(cmake/modules/CodeCoverage.cmake)
  SETUP_TARGET_FOR_COVERAGE(coverage ctest coverage)
endif(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# MFEM
if (NOT TARGET mfem)
  message(STATUS "Searching for MFEM")
  find_package(MFEM REQUIRED)
else(NOT TARGET mfem)
  message(STATUS "MFEM target already set")
endif(NOT TARGET mfem)

if (MFEM_USE_MPI)
  message(STATUS "MFEM: using MPI, so mfem_mgis will use MPI")
  find_package(MPI REQUIRED)
endif(MFEM_USE_MPI)

# MFontGenericInterface
find_package (MFrontGenericInterface REQUIRED)

# Support for OpenMP
if(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP) 
  find_package(OpenMP REQUIRED)
endif(MFEM_USE_OPENMP OR MFEM_USE_LEGACY_OPENMP)
 
# find mfront for testing purposes
include(cmake/modules/tfel.cmake)  

# documentation
if (NOT TARGET doc)
  add_custom_target(doc)
endif(NOT TARGET doc)

option(enable-website "enable generation of the website" ON)
if(enable-website)
  set(MFEMMGIS_GENERATE_WEBSITE ON)
  # Looking for pandoc (http://johnmacfarlane.net/pandoc)
  include(cmake/modules/pandoc.cmake)
  if(MFEMMGIS_HAVE_PANDOC)
    set(MFEMMGIS_GENERATE_WEBSITE ON)
  endif(MFEMMGIS_HAVE_PANDOC)
  if(MFEMMGIS_GENERATE_WEBSITE)
    list(APPEND CPACK_COMPONENTS_ALL website)
    set(CPACK_COMPONENT_WEBSITE_DESCRIPTION
      "Contains the MFEM/MGIS web site")
  endif(MFEMMGIS_GENERATE_WEBSITE)
else(enable-website)
  set(MFEMMGIS_GENERATE_WEBSITE OFF)
  message(STATUS "Generation of the MFEM/MGIS website disabled")
endif(enable-website)

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")

if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)

mfem_mgis_buildenv()

add_subdirectory(docs)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)
