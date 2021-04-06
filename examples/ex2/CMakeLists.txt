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
file(COPY ${MFEMMGIS_TESTDIR}/Plasticity.mfront
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_custom_target(mfront
  COMMAND mfront --obuild --interface=generic Plasticity.mfront
  COMMENT "Performing mfront code generation")
add_dependencies(Ssna303 mfront)


Include(ExternalProject)
ExternalProject_Add(ssna303_mesh_download
    PREFIX mesh
    URL https://icecube-eu-290.icedrive.io/download?p=_QX5pCKPMSnjByW6RnzM2yTxA8xs1QC4c2B.jUCE_u6lhrLk9nbi2fYCrTU2NNC_Dgg0pbW.Y2iAFsZt2WdVE3F9Zde1IFREqF.SUs9PStV4MjFF0a1eSAUHDFXwi0mCxNkJ84d4gbK7vdAAIRByVTvO5FpkoqpgVNtyMzkGniuuxm4xdKaxjOlhOT92sFcUb4aMhvrdvvauOL7vnfe5hA--
    URL_HASH SHA512=11a19d29e440c804930e6bf68ac81ce9233d40972e816daead0a27fefc0bc7e57c36f55578b5a123ccbdf15efb76df834536d9e9bf07b39704c0f19e4f41a0d5
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DOWNLOAD_NAME ssna303.msh
    DOWNLOAD_NO_EXTRACT 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    )
add_dependencies(Ssna303 ssna303_mesh_download)




