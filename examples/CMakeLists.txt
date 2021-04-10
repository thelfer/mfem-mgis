
set(UTILITY_FILE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UnitTestingUtilities.hxx)
add_custom_command(
  OUTPUT  ${UTILITY_FILE}
  COMMAND ${PROJECT_SOURCE_DIR}/examples/ex1/sed_ex1.sh ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx ${UTILITY_FILE}
  DEPENDS ${PROJECT_SOURCE_DIR}/examples/ex1/sed_ex1.sh ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx
  COMMENT "seding the UnitTestingUtilities.hxx")

add_custom_target(generate-ex1 ALL
  DEPENDS ${UTILITY_FILE})

install(FILES
  ${UTILITY_FILE}
  ${PROJECT_SOURCE_DIR}/tests/UniaxialTensileTest.cxx
  ex1/CMakeLists.txt
  ex1/Makefile
  DESTINATION share/mfem-mgis/examples/ex1)
