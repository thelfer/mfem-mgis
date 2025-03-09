install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
  DESTINATION share/doc/mfem-mgis
  COMPONENT doxygen)
install(FILES
  ${PROJECT_SOURCE_DIR}/docs/doxygen/img/background.svg
  DESTINATION share/doc/mfem-mgis/html/img/
  COMPONENT doxygen)
