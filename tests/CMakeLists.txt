function(add_test_case name)
  add_executable(${name} EXCLUDE_FROM_ALL ${name}.cxx)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${name} PRIVATE MFEMMGIS)
  add_test(NAME ${name} COMMAND ${name})
  add_dependencies(check ${name})
endfunction()

function(add_test_case_mpi name parameters)
  add_executable(${name} EXCLUDE_FROM_ALL ${name}.cxx)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${name} PRIVATE MFEMMGIS)
  add_test(NAME ${name} COMMAND ${name} ${${parameters}} "--parallel" "0")
  add_test(NAME ${name}P COMMAND mpirun -n 2 ./${name} ${${parameters}} "--parallel" "1")
  add_dependencies(check ${name})
endfunction()

add_test_case(UnitTestParameters)

if(MFEM_MGIS_HAVE_TFEL)
  mfront_behaviours_check_library(BehaviourTest
    Elasticity
    OrthotropicElasticity
    Plasticity
    Mazars
    SaintVenantKirchhoffElasticity
    SaintVenantKirchhoffElasticityV2
    MiehePhaseFieldDamage
    MicromorphicDamageI_SpectralSplit
    MicromorphicDamageII
    StationaryNonLinearHeatTransfer
    ImplicitNortonThreshold
    ElasticitySatohTest)

  set(SatohTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh" "--library" "$<TARGET_FILE:BehaviourTest>")
  add_test_case_mpi(SatohTest SatohTestParameters)

  add_executable(PeriodicTest
    EXCLUDE_FROM_ALL
    PeriodicTest.cxx)
  target_link_libraries(PeriodicTest
    PRIVATE MFEMMGIS)
  add_dependencies(check PeriodicTest)
  
  function(add_periodic_test ncase nsolver)
    set(test "PeriodicTest-${ncase}-${nsolver}")
    set(PeriodicTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh" "--library" "$<TARGET_FILE:BehaviourTest>" "--test-case" "${ncase}" "--linearsolver" "${nsolver}")
    if("1" STREQUAL ${nsolver} OR "0" STREQUAL ${nsolver})
    	add_test(NAME ${test}  COMMAND ./PeriodicTest ${PeriodicTestParameters} "--parallel" "0")
		endif()
    add_test(NAME ${test}P  COMMAND mpirun -n 2 ./PeriodicTest ${PeriodicTestParameters} "--parallel" "1")
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
        set_property(TEST ${test}P
          PROPERTY DEPENDS BehaviourTest
          PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
        set_property(TEST ${test}P
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_periodic_test)
  
  foreach(ncase RANGE 0 5)
    foreach(nsolver RANGE 0 5)
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
    set(UniaxialTensileTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh" "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "${behaviour}" "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref" "--internal-state-variable" "${internal_state_variable}")
    add_test(NAME ${test}  COMMAND ./UniaxialTensileTest ${UniaxialTensileTestParameters} "--parallel" "0")
    add_test(NAME ${test}P  COMMAND mpirun -n 2 ./UniaxialTensileTest ${UniaxialTensileTestParameters} "--parallel" "1")
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

  add_executable(ImposedPressureTest
    EXCLUDE_FROM_ALL
    ImposedPressureTest.cxx)
  target_include_directories(ImposedPressureTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(ImposedPressureTest
    PRIVATE MFEMMGIS)
  add_dependencies(check ImposedPressureTest)

  function(add_imposed_pressure_test behaviour)
    set(test "ImposedPressureTest-${behaviour}")
    set(ImposedPressureTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh" "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "${behaviour}" "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref")
    add_test(NAME ${test}  COMMAND ./ImposedPressureTest ${ImposedPressureTestParameters} "--parallel" "0")
    add_test(NAME ${test}P  COMMAND mpirun -n 2 ./ImposedPressureTest ${ImposedPressureTestParameters} "--parallel" "1")
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_imposed_pressure_test)

  add_imposed_pressure_test(Elasticity)

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
     COMMAND ./PartialQuadratureFunctionTest
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "${behaviour}"
     "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
     "--internal-state-variable" "${internal_state_variable}" "--parallel" 0)
    add_test(NAME ${test}P
     COMMAND mpirun -n 2 ./PartialQuadratureFunctionTest
     "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
     "--library" "$<TARGET_FILE:BehaviourTest>"
     "--behaviour" "${behaviour}"
     "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
     "--internal-state-variable" "${internal_state_variable}" "--parallel" 0)
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
    set(StationaryNonLinearHeatTransferTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"  "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "StationaryNonLinearHeatTransfer" "--linearsolver" "${nsolver}")
    add_test(NAME ${test} COMMAND ./StationaryNonLinearHeatTransferTest ${StationaryNonLinearHeatTransferTestParameters} "--parallel" "0")
    add_test(NAME ${test}P COMMAND mpirun -n 2 ./StationaryNonLinearHeatTransferTest ${StationaryNonLinearHeatTransferTestParameters} "--parallel" "1")
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

  add_executable(L2ProjectionTest
    EXCLUDE_FROM_ALL
    L2ProjectionTest.cxx)
  target_include_directories(L2ProjectionTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(L2ProjectionTest
    PRIVATE MFEMMGIS)
  add_dependencies(check L2ProjectionTest)
  
  function(add_l2_projection_test nsolver)
    set(test "L2ProjectionTest-${nsolver}")
    set(L2ProjectionTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver" "${nsolver}")
    add_test(NAME ${test}  COMMAND ./L2ProjectionTest ${L2ProjectionTestParameters} "--parallel" "0")
    add_test(NAME ${test}P  COMMAND mpirun -n 2 ./L2ProjectionTest ${L2ProjectionTestParameters} "--parallel" "1")
  endfunction(add_l2_projection_test)

  add_l2_projection_test(0)

  if(MFEM_USE_SUITESPARSE)
          message(STATUS "micromorphic_damage_2d has been disabled with SUITE SPARSE")
          #add_micromorphic_damage_2d_test(2)
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
    set(MicromorphicDamage2DTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "MiehePhaseFieldDamage" "--linearsolver" "${nsolver}")
    add_test(NAME ${test}  COMMAND ./MicromorphicDamage2DTest ${MicromorphicDamage2DTestParameters} "--parallel" "0")
    add_test(NAME ${test}P  COMMAND mpirun -n 2 ./MicromorphicDamage2DTest ${MicromorphicDamage2DTestParameters} "--parallel" "1")
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
          message(STATUS "micromorphic_damage_2d has been disabled with SUITE SPARSE")
          #add_micromorphic_damage_2d_test(2)
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
    set(MicromorphicDamage2DTest2Parameteres "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "MicromorphicDamageII" "--linearsolver" "${nsolver}")
    add_test(NAME ${test} COMMAND ./MicromorphicDamage2DTest2 ${MicromorphicDamage2DTest2Parameteres} "--parallel" "0")
    add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage2DTest2 ${MicromorphicDamage2DTest2Parameteres} "--parallel" "1")
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
    set(${test}Parameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar-with-defect-3D.msh" "--library" "$<TARGET_FILE:BehaviourTest>" "--behaviour" "MicromorphicDamageII" "--linearsolver" "${nsolver}")
    add_test(NAME ${test} COMMAND ./MicromorphicDamage3DTest ${${test}Parameters} "--parallel" "0")
    add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage3DTest ${${test}Parameters} "--parallel" "1")

    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
      set_property(TEST ${test}P
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}
        PROPERTY DEPENDS BehaviourTest)
      set_property(TEST ${test}P
        PROPERTY DEPENDS BehaviourTest)
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
  endfunction(add_micromorphic_damage_3d_test)
  
	add_micromorphic_damage_2d_test2(0)
	add_micromorphic_damage_2d_test(0)
#	add_micromorphic_damage_2d_test2(1)
	add_micromorphic_damage_3d_test(0)
#	add_micromorphic_damage_3d_test(1)
  if(MFEM_USE_SUITESPARSE)
#	add_micromorphic_damage_2d_test2(2)
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
  
#    add_test(NAME Restart_Parallel_Read_Mode_MPI_2
#      COMMAND mpirun -n 2 ./ParallelReadMode --mesh ${CMAKE_CURRENT_SOURCE_DIR}/mesh-test. --library "$<TARGET_FILE:BehaviourTest>" --read-mesh-mode Restart)
    
  endif(MFEM_USE_MPI)

  add_executable(PostProSubMesh
      EXCLUDE_FROM_ALL
      PostProSubMesh.cxx)
  target_link_libraries(PostProSubMesh
    PRIVATE MFEMMGIS)
  add_dependencies(check PostProSubMesh)
  
  if (MFEM_USE_MPI)
    add_test(NAME PostProSubMeshPar
      COMMAND mpirun -n 2 ./PostProSubMesh --mesh ${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh --library "$<TARGET_FILE:BehaviourTest>")
  else(MFEM_USE_MPI)
  add_test(NAME PostProSubMeshSeq
    COMMAND ./PostProSubMesh --mesh ${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh --library "$<TARGET_FILE:BehaviourTest>")
  endif(MFEM_USE_MPI)


## RVE examples

  if (MFEM_USE_MPI)
    set(RVEParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh" "--library" "$<TARGET_FILE:BehaviourTest>" "--post-processing" "1")
    add_test_case_mpi(RveNonLinearElastic RVEParameters)
    add_test_case_mpi(RveNonLinearElasticV2 RVEParameters)
    set(MixedOxideFuelsParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/inclusion.msh" "--order" "1" "--library" "$<TARGET_FILE:BehaviourTest>" "--post-processing" "1")
    add_test_case_mpi(MixedOxideFuels MixedOxideFuelsParameters)
  endif(MFEM_USE_MPI)

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
