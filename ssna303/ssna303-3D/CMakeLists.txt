cmake_minimum_required(VERSION 3.13)
# Define the executable to build for this example
add_executable(Ssna303 ssna303.cxx)
target_link_libraries(Ssna303 PRIVATE MFEMMGIS)

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/ssna303/ssna303-3D/Plasticity.mfront
  DESTINATION ${PROJECT_BINARY_DIR})

# Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Plasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Plasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-libBehaviour-Ssna303 ALL
  DEPENDS "src/libBehaviour.so")

add_dependencies(check Ssna303 generate-libBehaviour-Ssna303)

# Define a test that will be called during the "make check"
function(add_Ssna303_test_ex)
  set(test "Ssna303")
  add_test(NAME ${test} COMMAND Ssna303)
endfunction(add_Ssna303_test_ex)

add_Ssna303_test_ex()

