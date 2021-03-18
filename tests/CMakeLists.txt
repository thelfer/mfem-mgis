mfem_mgis_share(tests Plasticity.mfront)

mfront_behaviours_check_library(BehaviourTest
  Elasticity
  OrthotropicElasticity
  Plasticity
  Mazars
  SaintVenantKirchhoffElasticity)

add_executable(PeriodicTest
  EXCLUDE_FROM_ALL
  PeriodicTest.cxx)
target_link_libraries(PeriodicTest
  PRIVATE MFEMMGIS)
  add_dependencies(check PeriodicTest)

function(add_periodic_test ncase nsolver)
  set(test "PeriodicTest-${ncase}-${nsolver}")
  add_test(NAME ${test}
   COMMAND PeriodicTest 
   "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
   "--library" "$<TARGET_FILE:BehaviourTest>"
   "--test-case" "${ncase}"
   "--linearsolver" "${nsolver}")
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
endfunction(add_periodic_test)

foreach(ncase RANGE 0 5)
  foreach(nsolver RANGE 0 1)
    add_periodic_test(${ncase} ${nsolver})
  endforeach(nsolver)
endforeach(ncase)

add_executable(UniaxialTensileTest
  EXCLUDE_FROM_ALL
  UniaxialTensileTest.cxx)
target_link_libraries(UniaxialTensileTest
  PRIVATE MFEMMGIS)
add_dependencies(check UniaxialTensileTest)

function(add_unilateral_test behaviour internal_state_variable)
  set(test "UniaxialTensileTest-${behaviour}")
  add_test(NAME ${test}
   COMMAND UniaxialTensileTest
   "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
   "--library" "$<TARGET_FILE:BehaviourTest>"
   "--behaviour" "${behaviour}"
   "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
   "--internal-state-variable" "${internal_state_variable}")
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
endfunction(add_unilateral_test)

add_unilateral_test(OrthotropicElasticity EquivalentStrain)
add_unilateral_test(Plasticity EquivalentPlasticStrain)
add_unilateral_test(Mazars Damage)
add_unilateral_test(SaintVenantKirchhoffElasticity EquivalentStrain)

if (MFEM_USE_MPI)
  add_executable(PeriodicTestP
    EXCLUDE_FROM_ALL
    PeriodicTest.cxx)
  set_target_properties(PeriodicTestP PROPERTIES COMPILE_FLAGS "-DDO_USE_MPI" )
  target_link_libraries(PeriodicTestP
    PRIVATE MFEMMGIS)
  add_dependencies(check PeriodicTestP)
  
  function(add_periodic_testp ncase nsolver nbprocs)
    set(test "PeriodicTestP-${ncase}-${nsolver}-${nbprocs}")
    add_test(NAME ${test}
     COMMAND mpirun -n ${nbprocs} PeriodicTestP 
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--test-case" "${ncase}"
     "--linearsolver" "${nsolver}")
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_periodic_testp)
  
  foreach(ncase RANGE 0 0)
    foreach(nsolver RANGE 1 1)
      foreach(nbprocs RANGE 1 4 3)
	add_periodic_testp(${ncase} ${nsolver} ${nbprocs})
      endforeach(nbprocs)
    endforeach(nsolver)
  endforeach(ncase)

endif()
