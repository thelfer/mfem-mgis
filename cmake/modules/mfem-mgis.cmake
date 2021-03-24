set(EXPORT_INSTALL_PATH "share/mfem-mgis/cmake")

function(mfem_mgis_buildenv)
  set(METIS_DIR "$ENV{METIS_DIR}")
  set(HYPRE_DIR "$ENV{HYPRE_DIR}")
  if (MFEM_USE_MPI)
    set(MFEMMGIS_CXX "${MPI_CXX_COMPILER}")
  else()
    set(MFEMMGIS_CXX "${CMAKE_CXX_COMPILER}")
  endif()
  string(REGEX REPLACE "/mfront$" "" MFRONT_PATH "${MFRONT}")
  string(REGEX REPLACE "/[^/]*$" "" MPICXX_PATH "${MPI_CXX_COMPILER}")
  execute_process(COMMAND "which" "cmake" OUTPUT_VARIABLE CMAKE_CP_PATH)
  if(CMAKE_CP_PATH STREQUAL "")
    message( FATAL_ERROR  "Program cmake not found")
  endif()
  string(REGEX REPLACE "/cmake\n$" "" CMAKE_CP_PATH "${CMAKE_CP_PATH}")
  configure_file(cmake/env.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/env.sh  @ONLY)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/env.sh
    DESTINATION share/examples/)
  add_custom_target(sed_ex1 ALL
    COMMAND  ${CMAKE_SOURCE_DIR}/examples/ex1/sed_ex1.sh ${CMAKE_SOURCE_DIR}/tests/UniaxialTensileTest.cxx  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UniaxialTensileTest.cxx
    COMMENT "seding the UniaxialTensileTest.cxx"
    DEPENDS MFEMMGIS)
endfunction(mfem_mgis_buildenv)

function(mfem_mgis_header dir file)
  install(FILES ${dir}/${file}
    DESTINATION "include/${dir}")
endfunction(mfem_mgis_header)

function(mfem_mgis_share dir file)
  install(FILES ${file}
    DESTINATION "share/${dir}")
endfunction(mfem_mgis_share)

function(mfem_mgis_library name)
  if(${ARGC} LESS 2)
    message(FATAL_ERROR "mfem_mgis_library_internal : no source specified")
  endif(${ARGC} LESS 2)
  add_library(${name} STATIC ${ARGN})
  target_include_directories(${name}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PUBLIC $<INSTALL_INTERFACE:include>
  )
  target_include_directories(${name}
    SYSTEM
    PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
    PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
  target_link_libraries(${name} 
    PUBLIC mgis::MFrontGenericInterface
    PUBLIC mfem)
  if(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION bin)
  else(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION lib${LIB_SUFFIX})
  endif(WIN32)
  install(EXPORT ${name} DESTINATION ${EXPORT_INSTALL_PATH}
    NAMESPACE mfem-mgis:: FILE ${name}Config.cmake)
  
endfunction(mfem_mgis_library)
