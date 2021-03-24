
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UniaxialTensileTest.cxx
  DESTINATION share/examples/ex1)

install(FILES
  ex1/CMakeLists.txt
  ex1/Makefile
  DESTINATION share/examples/ex1)
