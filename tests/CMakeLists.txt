mfront_behaviours_check_library(BehaviourTest
  Elasticity)

add_executable(PeriodicTest
  EXCLUDE_FROM_ALL
  PeriodicTest.cxx)
target_link_libraries(PeriodicTest
  PRIVATE MFEMMGIS)

function(add_periodic_test ncase)
  set(test "PeriodicTest-${ncase}")
  add_test(NAME ${test}
   COMMAND PeriodicTest 
   "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
   "--library" "$<TARGET_FILE:BehaviourTest>"
   "--test-case" "${ncase}")
  add_dependencies(check PeriodicTest)
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
endfunction(add_periodic_test)

add_periodic_test(0)
add_periodic_test(1)
add_periodic_test(2)
add_periodic_test(3)
add_periodic_test(4)
add_periodic_test(5)

# add_executable(BehaviourIntegratorTest
#   EXCLUDE_FROM_ALL
#   BehaviourIntegratorTest.cxx)
# target_link_libraries(BehaviourIntegratorTest
#   PRIVATE MFEMMGIS)
