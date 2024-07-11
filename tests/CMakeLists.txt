if(MFEM_MGIS_HAVE_TFEL)
  mfront_behaviours_check_library(BehaviourTest
    Elasticity
    OrthotropicElasticity
    Plasticity
    Mazars
    SaintVenantKirchhoffElasticity
    MiehePhaseFieldDamage
    MicromorphicDamageI_SpectralSplit
    MicromorphicDamageII
    StationaryNonLinearHeatTransfer)
  
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
  target_include_directories(UniaxialTensileTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(UniaxialTensileTest
    PRIVATE MFEMMGIS)
  add_dependencies(check UniaxialTensileTest)
  
  function(add_uniaxial_tensile_test behaviour internal_state_variable)
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
  endfunction(add_uniaxial_tensile_test)
  
  add_uniaxial_tensile_test(OrthotropicElasticity EquivalentStrain)
  add_uniaxial_tensile_test(Plasticity EquivalentPlasticStrain)
  add_uniaxial_tensile_test(Mazars Damage)
  add_uniaxial_tensile_test(SaintVenantKirchhoffElasticity EquivalentStrain)

  add_executable(PartialQuadratureFunctionTest
    EXCLUDE_FROM_ALL
    PartialQuadratureFunctionTest.cxx)
  target_include_directories(PartialQuadratureFunctionTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(PartialQuadratureFunctionTest
    PRIVATE MFEMMGIS)
  add_dependencies(check PartialQuadratureFunctionTest)
  
  function(add_partial_quadrature_function_test behaviour internal_state_variable)
    set(test "PartialQuadratureFunctionTest-${behaviour}")
    add_test(NAME ${test}
     COMMAND PartialQuadratureFunctionTest
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
  endfunction(add_partial_quadrature_function_test)
  
  add_partial_quadrature_function_test(OrthotropicElasticity EquivalentStrain)
  
  add_executable(StationaryNonLinearHeatTransferTest
    EXCLUDE_FROM_ALL
    StationaryNonLinearHeatTransferTest.cxx)
  target_include_directories(StationaryNonLinearHeatTransferTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(StationaryNonLinearHeatTransferTest
    PRIVATE MFEMMGIS)
  add_dependencies(check StationaryNonLinearHeatTransferTest)
  
  function(add_stationary_non_linear_heat_transfer_test nsolver)
    set(test "StationaryNonLinearHeatTransferTest-${nsolver}")
    add_test(NAME ${test}
     COMMAND StationaryNonLinearHeatTransferTest
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "StationaryNonLinearHeatTransfer"
     "--linearsolver" "${nsolver}"
#     "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
	)
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_stationary_non_linear_heat_transfer_test)

  add_stationary_non_linear_heat_transfer_test(1)
  if(MFEM_USE_SUITESPARSE)
	add_stationary_non_linear_heat_transfer_test(2)
  endif(MFEM_USE_SUITESPARSE)

  add_executable(MicromorphicDamage2DTest
    EXCLUDE_FROM_ALL
    MicromorphicDamage2DTest.cxx)
  target_include_directories(MicromorphicDamage2DTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(MicromorphicDamage2DTest
    PRIVATE MFEMMGIS)
  add_dependencies(check MicromorphicDamage2DTest)
  
  function(add_micromorphic_damage_2d_test nsolver)
    set(test "MicromorphicDamage2DTest-${nsolver}")
    add_test(NAME ${test}
     COMMAND ./MicromorphicDamage2DTest
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "MiehePhaseFieldDamage"
     "--linearsolver" "${nsolver}"
#     "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
	)
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_micromorphic_damage_2d_test)

  if(MFEM_USE_SUITESPARSE)
	add_micromorphic_damage_2d_test(2)
  endif(MFEM_USE_SUITESPARSE)

  add_executable(MicromorphicDamage2DTest2
    EXCLUDE_FROM_ALL
    MicromorphicDamage2DTest2.cxx)
  target_include_directories(MicromorphicDamage2DTest2
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(MicromorphicDamage2DTest2
    PRIVATE MFEMMGIS)
  add_dependencies(check MicromorphicDamage2DTest2)

  function(add_micromorphic_damage_2d_test2 nsolver)
    set(test "MicromorphicDamage2DTest2-${nsolver}")
    add_test(NAME ${test}
     COMMAND ./MicromorphicDamage2DTest2
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "MicromorphicDamageII"
     "--linearsolver" "${nsolver}"
#     "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
	)
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_micromorphic_damage_2d_test2)

  add_executable(MicromorphicDamage3DTest
    EXCLUDE_FROM_ALL
    MicromorphicDamage3DTest.cxx)
  target_include_directories(MicromorphicDamage3DTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(MicromorphicDamage3DTest
    PRIVATE MFEMMGIS)
  add_dependencies(check MicromorphicDamage3DTest)

  function(add_micromorphic_damage_3d_test nsolver)
    set(test "MicromorphicDamage3DTest-${nsolver}")
    add_test(NAME ${test}
     COMMAND MicromorphicDamage3DTest
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar-with-defect-3D.msh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "MicromorphicDamageII"
     "--linearsolver" "${nsolver}"
	)
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_micromorphic_damage_3d_test)


  
  if(MFEM_USE_SUITESPARSE)
	add_micromorphic_damage_2d_test2(2)
#	add_micromorphic_damage_3d_test(2)
  endif(MFEM_USE_SUITESPARSE)

  if (MFEM_USE_MPI)

    add_executable(ParallelReadMode
      EXCLUDE_FROM_ALL
      ParallelReadMode.cxx)
    set_target_properties(ParallelReadMode PROPERTIES COMPILE_FLAGS "-DDO_USE_MPI" )
    target_link_libraries(ParallelReadMode
      PRIVATE MFEMMGIS)
    add_dependencies(check ParallelReadMode)
  
    add_test(NAME Restart_Parallel_Read_Mode_MPI_4
      COMMAND mpirun -n 4 ./ParallelReadMode --mesh ${CMAKE_CURRENT_SOURCE_DIR}/mesh-explorer.mesh. --library "$<TARGET_FILE:BehaviourTest>" --read-mesh-mode Restart)
      

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
       COMMAND mpirun -n ${nbprocs} ./PeriodicTestP 
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

    if(MFEM_USE_MUMPS)

      add_periodic_testp(0 3 1)
      add_periodic_testp(0 3 4)

      add_executable(UniaxialTensileTestP
      	EXCLUDE_FROM_ALL
      	UniaxialTensileTest.cxx)
      set_target_properties(UniaxialTensileTestP PROPERTIES COMPILE_FLAGS "-DDO_USE_MPI" )
      target_include_directories(UniaxialTensileTestP
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
      target_link_libraries(UniaxialTensileTestP
	    PRIVATE MFEMMGIS)
      add_dependencies(check UniaxialTensileTestP)
      
      function(add_uniaxial_tensile_testp behaviour internal_state_variable linear_solver nbprocs)
	    set(test "UniaxialTensileTestP-${behaviour}-${linear_solver}-${nbprocs}")
	    add_test(NAME ${test}
	      COMMAND mpirun -n ${nbprocs} UniaxialTensileTestP
	      "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
	      "--library" "$<TARGET_FILE:BehaviourTest>"
	      "--linearsolver" "${linear_solver}"
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
      endfunction(add_uniaxial_tensile_testp)
      
      add_uniaxial_tensile_testp(OrthotropicElasticity EquivalentStrain 3 1)
      add_uniaxial_tensile_testp(OrthotropicElasticity EquivalentStrain 3 4)

    endif(MFEM_USE_MUMPS)
    
  endif(MFEM_USE_MPI)
else(MFEM_MGIS_HAVE_TFEL)
endif(MFEM_MGIS_HAVE_TFEL)

# Configure and export example 1
# This example is based on UniaxialTensileTest from the unit tests

set(EX1_DIR "share/mfem-mgis/examples/ex1")
install(FILES Plasticity.mfront cube.mesh 
  DESTINATION ${EX1_DIR})

install(FILES ex1CMakeLists.txt
  DESTINATION ${EX1_DIR}
  RENAME CMakeLists.txt)

install(FILES ex1Makefile
  DESTINATION ${EX1_DIR}
  RENAME Makefile)

set(UTILITY_FILE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UnitTestingUtilities.hxx)
add_custom_command(
  OUTPUT  ${UTILITY_FILE}
  COMMAND ${PROJECT_SOURCE_DIR}/tests/sed_utility.sh ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx ${UTILITY_FILE}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests/sed_utility.sh ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx
  COMMENT "seding the UnitTestingUtilities.hxx")

add_custom_target(generate-ex1 ALL
  DEPENDS ${UTILITY_FILE})

install(FILES
  ${UTILITY_FILE}
  ${PROJECT_SOURCE_DIR}/tests/UniaxialTensileTest.cxx
  DESTINATION ${EX1_DIR})

set(EX2_DIR "share/mfem-mgis/examples/ex2")
install(FILES PeriodicTest.cxx Elasticity.mfront cube_2mat_per.mesh
  DESTINATION ${EX2_DIR})

install(FILES ex2CMakeLists.txt
  DESTINATION ${EX2_DIR}
  RENAME CMakeLists.txt)
