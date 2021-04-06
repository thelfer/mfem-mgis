cmake_minimum_required(VERSION 3.13)
project(UnilateralTensileTestExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MFEM REQUIRED)
find_package(MFrontGenericInterface REQUIRED)
find_package(MFEMMGIS REQUIRED)

add_executable(Ssna303 ssna303.cxx)
target_link_libraries(Ssna303 mfem-mgis::MFEMMGIS)
set(MFEMMGIS_TESTDIR "${MFEMMGIS_DIR}/../examples/ex2")
file(COPY
  ${MFEMMGIS_TESTDIR}/Plasticity.mfront
  ${MFEMMGIS_TESTDIR}/ssna303.msh
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_custom_target(mfront
  COMMAND mfront --obuild --interface=generic Plasticity.mfront
  COMMENT "Performing mfront code generation")
add_dependencies(Ssna303 mfront)





