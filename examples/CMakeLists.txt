#  const char* mesh_file = "cube.mesh";
#  const char* behaviour = "Plasticity";
#  const char* library = "src/libBehaviour.so";
#  const char* reference_file = "Plasticity.ref";
#  const char* isv_name = "EquivalentPlasticStrain";

add_custom_target(sed_ex1 ALL
  COMMAND  ${CMAKE_SOURCE_DIR}/examples/ex1/sed_ex1.sh ${CMAKE_SOURCE_DIR}/tests/UniaxialTensileTest.cxx  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UniaxialTensileTest.cxx
  COMMENT "seding the UniaxialTensileTest.cxx")

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UniaxialTensileTest.cxx
  DESTINATION share/examples/ex1)

install(FILES
  ex1/CMakeLists.txt
  DESTINATION share/examples/ex1)
