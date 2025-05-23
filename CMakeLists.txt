cmake_minimum_required(VERSION 3.8)
project(mfem-mgis
  VERSION 0.2
  DESCRIPTION "The aim of the MFEM/MGIS project is to provide a C++ library to build advanced mechanical simulation based on MFront behaviours with MFEM.")

include(cmake/modules/mfem-mgis.cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# (must be placed *before* any add_subdirectory, cmake bug ?)
include(CTest)
enable_testing()

# portable-build
option(enable-portable-build "Produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)
option(BENCHMARKS "Compile benchmarks" OFF)

#compiler options
include(cmake/modules/compiler.cmake)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(cmake/modules/CodeCoverage.cmake)
  SETUP_TARGET_FOR_COVERAGE(coverage ctest coverage)
endif(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# MFEM
find_package(MFEM REQUIRED)

if (MFEM_USE_MPI)
  message(STATUS "MFEM: using MPI, so mfem_mgis will use MPI")
  find_package(MPI REQUIRED)
  # The following include directories are automatically filled within FindMFEM.cmake
  include_directories(SYSTEM ${HYPRE_INCLUDE_DIRS})
  include_directories(SYSTEM ${METIS_INCLUDE_DIRS})
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

# add a target to generate API documentation with Doxygen
option(enable-doxygen-doc "enable generation of the Doxygen documentation" OFF)
if(enable-doxygen-doc)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    add_custom_target(doc-html
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM)
    add_dependencies(doc doc-html)
    list(APPEND CPACK_COMPONENTS_ALL doxygen)
    set(CPACK_COMPONENT_DOXYGEN_DESCRIPTION
      "Contains the TFEL doxygen documentation")
    set(GENERATE_DOXYGEN ON)
  else(DOXYGEN_FOUND)
    message(FATAL_ERROR "doxygen is required")
  endif(DOXYGEN_FOUND)
else(enable-doxygen-doc)
  message(STATUS "doxygen documentation is not enabled")
endif(enable-doxygen-doc)

add_subdirectory(docs)

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")

if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test --output-on-failure -C $<CONFIGURATION>)
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test --output-on-failure -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test --output-on-failure )
   else(NOT TARGET check)
      add_custom_target(mgis_check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test --output-on-failure )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)

# HYPRE version
include(cmake/modules/hypre.cmake)  
set(MFEM_HYPRE_VERSION ${HYPRE_VERSION})
add_compile_definitions(MFEM_HYPRE_VERSION=${MFEM_HYPRE_VERSION})


mfem_mgis_buildenv()

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)
if(BENCHMARKS)
  add_subdirectory(benchmarks)
endif(BENCHMARKS)
