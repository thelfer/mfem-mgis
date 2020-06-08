cmake_minimum_required(VERSION 3.7)
project(mfem-mgis)

include(cmake/modules/mfem-mgis.cmake)

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

#compiler options
include(cmake/modules/compiler.cmake)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(cmake/modules/CodeCoverage.cmake)
  SETUP_TARGET_FOR_COVERAGE(coverage ctest coverage)
endif(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# MFEM
find_package(MFEM REQUIRED)

# MFontGenericInterface
find_package (MFrontGenericInterface REQUIRED)
 
# find mfront for testing purposes
include(cmake/modules/tfel.cmake)  

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)