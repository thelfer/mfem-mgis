mfront_behaviours_check_library(BehaviourTest
  Elasticity)

add_executable(PeriodicTest
  EXCLUDE_FROM_ALL
  PeriodicTest.cxx)
target_link_libraries(PeriodicTest
  PRIVATE MFEMMGIS)

add_test(NAME PeriodicTest
 COMMAND PeriodicTest 
 "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/square_2mat_per.mesh"
 "--library" "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check PeriodicTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST PeriodicTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST PeriodicTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

# add_executable(BehaviourIntegratorTest
#   EXCLUDE_FROM_ALL
#   BehaviourIntegratorTest.cxx)
# target_link_libraries(BehaviourIntegratorTest
#   PRIVATE MFEMMGIS)
