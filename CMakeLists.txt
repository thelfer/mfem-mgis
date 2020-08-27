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

# documentation
add_custom_target(doc)

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

add_subdirectory(docs)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)
