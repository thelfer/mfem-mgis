cmake_minimum_required(VERSION 3.13)
project(UnilateralTensileTestExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MFEM REQUIRED)
find_package(MFrontGenericInterface REQUIRED)
find_package(MFEMMGIS REQUIRED)

add_executable(ssna303 ssna303.cxx)
target_link_libraries(ssna303 mfem-mgis::MFEMMGIS)





