include(CMakePackageConfigHelpers)
set(EXPORT_INSTALL_PATH "share/mfem-mgis/cmake")

function(mfem_mgis_buildenv)
  if (MFEM_USE_MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
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
  get_target_property(MFEM_LIB_LIST MFEM::mfem INTERFACE_LINK_LIBRARIES)
  get_target_property(MFEM_INCLUDE_DIRS MFEM::mfem INTERFACE_INCLUDE_DIRECTORIES)
  configure_file(cmake/env.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/env.sh  @ONLY)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/env.sh
    DESTINATION share/mfem-mgis/examples/)
endfunction(mfem_mgis_buildenv)

function(mfem_mgis_header dir file)
  install(FILES ${dir}/${file}
    DESTINATION "include/${dir}")
endfunction(mfem_mgis_header)

function(mfem_mgis_library name)
  if(${ARGC} LESS 2)
    message(FATAL_ERROR "mfem_mgis_library_internal : no source specified")
  endif(${ARGC} LESS 2)
  add_library(${name} STATIC ${ARGN})
  target_include_directories(${name}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PUBLIC $<INSTALL_INTERFACE:include>
  )
  target_link_libraries(${name} 
    PUBLIC mgis::MFrontGenericInterface
    PUBLIC MFEM::mfem)
  if(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION bin)
  else(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION lib${LIB_SUFFIX})
  endif(WIN32)
  install(EXPORT ${name} DESTINATION ${EXPORT_INSTALL_PATH}
          EXPORT_LINK_INTERFACE_LIBRARIES
          NAMESPACE mfem-mgis:: FILE ${name}Targets.cmake)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name}Config.cmake.in)
    set(_package_config_file ${CMAKE_CURRENT_SOURCE_DIR}/${name}Config.cmake.in)
  else()
    set(_package_config_file ${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake.in)
    file(WRITE ${_package_config_file}
         "@PACKAGE_INIT@\n"
         "\n"
         "include(\"\${CMAKE_CURRENT_LIST_DIR}/${name}Targets.cmake\")")
  endif()
  # generate the config file that includes the exports
  configure_package_config_file(${_package_config_file}
                                "${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake"
                                INSTALL_DESTINATION ${EXPORT_INSTALL_PATH}
                                NO_SET_AND_CHECK_MACRO
                                NO_CHECK_REQUIRED_COMPONENTS_MACRO)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake"
          DESTINATION ${EXPORT_INSTALL_PATH})
endfunction(mfem_mgis_library)
