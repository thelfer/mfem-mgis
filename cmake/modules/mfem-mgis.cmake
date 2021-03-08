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
  target_include_directories(${name}
    SYSTEM PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>")
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
  if(MFEMMGIS_APPEND_SUFFIX)
    set(export_install_path "share/mfem-mgis-${MFEMMGIS_SUFFIX}/cmake")
  else(MFEMMGIS_APPEND_SUFFIX)
    set(export_install_path "share/mfem-mgis/cmake")
  endif(MFEMMGIS_APPEND_SUFFIX)
  install(EXPORT ${name} DESTINATION ${export_install_path}
    NAMESPACE mfem-mgis:: FILE ${name}Config.cmake)
  
endfunction(mfem_mgis_library)
