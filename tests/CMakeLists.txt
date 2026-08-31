# Helper function to create a test executable and add it to build-tests
function(add_test_executable_base name libs)
  add_executable(${name} EXCLUDE_FROM_ALL ${name}.cxx)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${name} PRIVATE "${libs}")
  add_dependencies(build-tests ${name})
endfunction()

function(add_test_executable name)
  add_test_executable_base("${name}" "mfem-mgis::MFEMMGIS")
endfunction()

function(add_tfel_test_executable name)
  if(MFEM_MGIS_ENABLE_TFEL_TESTS)
    add_test_executable_base("${name}" "mfem-mgis::MFEMMGIS;tfel::TFELTests")
  endif(MFEM_MGIS_ENABLE_TFEL_TESTS)
endfunction()

# Helper function to add a simple test case (serial only)
function(add_simple_test name)
  add_test_executable(${name})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# Helper function to add a test case with MPI support Usage: add_mpi_test(name
# parameters [NO_PARALLEL])
function(add_mpi_test name parameters)
  if(MFEM_USE_MPI)
    add_test_executable(${name})
    add_test(NAME ${name} COMMAND ${name} ${${parameters}} "--parallel" "0")
    add_test(NAME ${name}P COMMAND mpirun -n 2 ./${name} ${${parameters}}
                                   "--parallel" "1")
  endif(MFEM_USE_MPI)
endfunction()

# Simple test case (no TFEL dependency)
add_simple_test(UnitTestParameters)

if(MFEM_MGIS_HAVE_TFEL)
  # Compile MFront behaviours library
  mfront_behaviours_check_library(
    BehaviourTest
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
    ElasticitySatohTest
    UO2Shrinkage_RAPHAEL2008)

  # SatohTest with MPI support
  set(SatohTestParameters "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube2d.msh"
                          "--library" "$<TARGET_FILE:BehaviourTest>")
  add_mpi_test(SatohTest SatohTestParameters)

  # PeriodicTest executable
  add_test_executable(PeriodicTest)

  # Helper function for PeriodicTest cases
  function(add_periodic_test ncase nsolver)
    set(test "PeriodicTest-${ncase}-${nsolver}")
    set(PeriodicTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--test-case"
        "${ncase}"
        "--linearsolver"
        "${nsolver}")
    if("1" STREQUAL ${nsolver} OR "0" STREQUAL ${nsolver})
      add_test(NAME ${test} COMMAND ./PeriodicTest ${PeriodicTestParameters}
                                    "--parallel" "0")
      set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    endif()
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PeriodicTest
                                     ${PeriodicTestParameters} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # Generate all PeriodicTest cases (6x6 matrix)
  foreach(ncase RANGE 0 5)
    foreach(nsolver RANGE 0 5)
      add_periodic_test(${ncase} ${nsolver})
    endforeach()
  endforeach()

  # UniaxialTensileTest executable
  add_test_executable(UniaxialTensileTest)

  # Helper function for UniaxialTensileTest cases
  function(add_uniaxial_tensile_test behaviour internal_state_variable)
    set(test "UniaxialTensileTest-${behaviour}")
    set(UniaxialTensileTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test}
             COMMAND ./UniaxialTensileTest ${UniaxialTensileTestParameters}
                     "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./UniaxialTensileTest
                       ${UniaxialTensileTestParameters} "--parallel" "1")
    endif(MFEM_USE_MPI)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  endfunction()

  # UniaxialTensileTest cases
  add_uniaxial_tensile_test(OrthotropicElasticity EquivalentStrain)
  add_uniaxial_tensile_test(Plasticity EquivalentPlasticStrain)
  add_uniaxial_tensile_test(Mazars Damage)
  add_uniaxial_tensile_test(SaintVenantKirchhoffElasticity EquivalentStrain)

  # ImposedPressureTest executable
  add_test_executable(ImposedPressureTest)

  # Helper function for ImposedPressureTest cases
  function(add_imposed_pressure_test behaviour)
    set(test "ImposedPressureTest-${behaviour}")
    set(ImposedPressureTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref")
    add_test(NAME ${test}
             COMMAND ./ImposedPressureTest ${ImposedPressureTestParameters}
                     "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImposedPressureTest
                       ${ImposedPressureTestParameters} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # ImposedPressureTest2 executable
  add_test_executable(ImposedPressureTest2)

  # Helper function for ImposedPressureTest2 cases
  function(add_imposed_pressure_test2 behaviour)
    set(test "ImposedPressureTest2-${behaviour}")
    set(ImposedPressureTestParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref")
    add_test(NAME ${test}
             COMMAND ./ImposedPressureTest2 ${ImposedPressureTestParameters}
                     "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImposedPressureTest2
                       ${ImposedPressureTestParameters} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_imposed_pressure_test2(Elasticity)

  # PartialQuadratureFunctionTest executable
  add_tfel_test_executable(PartialQuadratureFunctionTest)

  # PartialQuadratureSpaceTest executable
  add_test_executable(PartialQuadratureSpaceTest)

  # Helper function for PartialQuadratureSpaceTest cases
  function(add_partial_quadrature_space_test behaviour internal_state_variable)
    set(test "PartialQuadratureSpaceTest-${behaviour}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test} COMMAND ./PartialQuadratureSpaceTest ${params}
                                  "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PartialQuadratureSpaceTest
                                     ${params} "--parallel" 1)
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_partial_quadrature_space_test(OrthotropicElasticity EquivalentStrain)

  # Helper function for PartialQuadratureFunctionTest cases
  function(add_partial_quadrature_function_test behaviour
           internal_state_variable)
    set(test "PartialQuadratureFunctionTest-${behaviour}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "${behaviour}"
        "--reference-file"
        "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
        "--internal-state-variable"
        "${internal_state_variable}")
    add_test(NAME ${test} COMMAND ./PartialQuadratureFunctionTest ${params}
                                  "--parallel" 0)
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./PartialQuadratureFunctionTest ${params}
                       "--parallel" 1)
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  add_partial_quadrature_function_test(OrthotropicElasticity EquivalentStrain)

  # StationaryNonLinearHeatTransferTest executable
  add_test_executable(StationaryNonLinearHeatTransferTest)

  # Helper function for StationaryNonLinearHeatTransferTest cases
  function(add_stationary_non_linear_heat_transfer_test nsolver)
    set(test "StationaryNonLinearHeatTransferTest-${nsolver}")
    set(params
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--behaviour"
        "StationaryNonLinearHeatTransfer"
        "--linearsolver"
        "${nsolver}")
    add_test(NAME ${test} COMMAND ./StationaryNonLinearHeatTransferTest
                                  ${params} "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./StationaryNonLinearHeatTransferTest
                       ${params} "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # StationaryNonLinearHeatTransferTest cases
  add_stationary_non_linear_heat_transfer_test(1)
  if(MFEM_USE_SUITESPARSE)
    add_stationary_non_linear_heat_transfer_test(2)
  endif()

  # PointWiseModelTest executable
  add_test_executable(PointWiseModelTest)

  # Helper function for PointWiseModelTest
  function(add_point_wise_model_test)
    set(test "PointWiseModelTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./PointWiseModelTest ${params} "--parallel"
                                  "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PointWiseModelTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  # PointWiseModelTest-2 executable
  add_test_executable(PointWiseModelTest-2)

  # Helper function for PointWiseModelTest-2
  function(add_point_wise_model_test_2)
    set(test "PointWiseModelTest-2")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./PointWiseModelTest-2 ${params} "--parallel"
                                  "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./PointWiseModelTest-2
                                     ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_point_wise_model_test_2()

  # L2ProjectionTest executable
  add_test_executable(L2ProjectionTest)

  # Helper function for L2ProjectionTest cases
  function(add_l2_projection_test nsolver)
    set(test "L2ProjectionTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./L2ProjectionTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./L2ProjectionTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_l2_projection_test(0)

  # GridFunctionTest executable
  add_test_executable(GridFunctionTest)

  # Helper function for GridFunctionTest
  function(add_grid_function_test)
    set(test "GridFunctionTest")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh")
    add_test(NAME ${test} COMMAND ./GridFunctionTest ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 ./GridFunctionTest ${params}
                                     "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_grid_function_test()

  # ImplicitGradientRegularizationTest executable
  add_test_executable(ImplicitGradientRegularizationTest)

  # Helper function for ImplicitGradientRegularizationTest cases
  function(add_implicit_gradient_regularization_test nsolver)
    set(test "ImplicitGradientRegularizationTest-${nsolver}")
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh" "--linearsolver"
               "${nsolver}")
    add_test(NAME ${test} COMMAND ./ImplicitGradientRegularizationTest
                                  ${params} "--parallel" "0")
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P
               COMMAND mpirun -n 2 ./ImplicitGradientRegularizationTest
                       ${params} "--parallel" "1")
    endif(MFEM_USE_MPI)
  endfunction()

  add_implicit_gradient_regularization_test(0)

  # # MicromorphicDamage2DTest executable
  # add_test_executable(MicromorphicDamage2DTest)

  # # Helper function for MicromorphicDamage2DTest cases
  # function(add_micromorphic_damage_2d_test nsolver)
  #   set(test "MicromorphicDamage2DTest-${nsolver}")
  #   set(params
  #       "--mesh"
  #       "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
  #       "--library"
  #       "$<TARGET_FILE:BehaviourTest>"
  #       "--behaviour"
  #       "MiehePhaseFieldDamage"
  #       "--linearsolver"
  #       "${nsolver}")
  #   add_test(NAME ${test} COMMAND ./MicromorphicDamage2DTest ${params}
  #                                 "--parallel" "0")
  #   set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  #   if(MFEM_USE_MPI)
  #     add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage2DTest
  #                                    ${params} "--parallel" "1")
  #     set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
  #   endif(MFEM_USE_MPI)
  # endfunction()

  # # MicromorphicDamage2DTest2 executable
  # add_test_executable(MicromorphicDamage2DTest2)

  # # Helper function for MicromorphicDamage2DTest2 cases
  # function(add_micromorphic_damage_2d_test2 nsolver)
  #   set(test "MicromorphicDamage2DTest2-${nsolver}")
  #   set(params
  #       "--mesh"
  #       "${CMAKE_CURRENT_SOURCE_DIR}/bar.msh"
  #       "--library"
  #       "$<TARGET_FILE:BehaviourTest>"
  #       "--behaviour"
  #       "MicromorphicDamageII"
  #       "--linearsolver"
  #       "${nsolver}")
  #   add_test(NAME ${test} COMMAND ./MicromorphicDamage2DTest2 ${params}
  #                                 "--parallel" "0")
  #   set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  #   if(MFEM_USE_MPI)
  #     add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage2DTest2
  #                                    ${params} "--parallel" "1")
  #     set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
  #   endif(MFEM_USE_MPI)
  # endfunction()

  # # MicromorphicDamage3DTest executable
  # add_test_executable(MicromorphicDamage3DTest)

  # # Helper function for MicromorphicDamage3DTest cases
  # function(add_micromorphic_damage_3d_test nsolver)
  #   set(test "MicromorphicDamage3DTest-${nsolver}")
  #   set(params
  #       "--mesh"
  #       "${CMAKE_CURRENT_SOURCE_DIR}/bar-with-defect-3D.msh"
  #       "--library"
  #       "$<TARGET_FILE:BehaviourTest>"
  #       "--behaviour"
  #       "MicromorphicDamageII"
  #       "--linearsolver"
  #       "${nsolver}")
  #   add_test(NAME ${test} COMMAND ./MicromorphicDamage3DTest ${params}
  #                                 "--parallel" "0")
  #   set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
  #   if(MFEM_USE_MPI)
  #     add_test(NAME ${test}P COMMAND mpirun -n 2 ./MicromorphicDamage3DTest
  #                                    ${params} "--parallel" "1")
  #     set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
  #   endif(MFEM_USE_MPI)
  # endfunction()

  # # Micromorphic damage tests
  # add_micromorphic_damage_2d_test2(0)
  # add_micromorphic_damage_2d_test(0)
  # add_micromorphic_damage_3d_test(0)

  # Helper function for PredictionTest cases
  function(add_prediction_test test)
    add_test_executable(${test})
    set(params "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh" "--library"
               "$<TARGET_FILE:BehaviourTest>")
    add_test(NAME ${test} COMMAND "./${test}" ${params} "--parallel" "0")
    set_property(TEST ${test} PROPERTY DEPENDS BehaviourTest)
    if(MFEM_USE_MPI)
      add_test(NAME ${test}P COMMAND mpirun -n 2 "./${test}" ${params}
                                     "--parallel" "1")
      set_property(TEST ${test}P PROPERTY DEPENDS BehaviourTest)
    endif(MFEM_USE_MPI)
  endfunction()

  # PredictionTest cases
  add_prediction_test(PredictionTest)
  add_prediction_test(PredictionTest2)

  # ParallelReadMode executable (MPI only)
  if(MFEM_USE_MPI)
    add_test_executable(ParallelReadMode)
    set_target_properties(ParallelReadMode PROPERTIES COMPILE_FLAGS
                                                      "-DDO_USE_MPI")
    # Commented out test - can be enabled if needed add_test(NAME
    # Restart_Parallel_Read_Mode_MPI_2 COMMAND mpirun -n 2 ./ParallelReadMode
    # --mesh ${CMAKE_CURRENT_SOURCE_DIR}/mesh-test. --library
    # "$<TARGET_FILE:BehaviourTest>" --read-mesh-mode Restart)
  endif()

  # PostProSubMesh executable (MPI only)
  if(MFEM_USE_MPI)
    add_test_executable(PostProSubMesh)
    add_test(
      NAME PostProSubMeshPar
      COMMAND
        mpirun -n 2 ./PostProSubMesh --mesh
        ${CMAKE_CURRENT_SOURCE_DIR}/beam-tet.mesh --library
        "$<TARGET_FILE:BehaviourTest>")
  endif()

  # RVE examples

  if(MFEM_USE_MPI)
    set(RVEParameters
        "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh" "--library"
        "$<TARGET_FILE:BehaviourTest>" "--post-processing" "1")
    add_mpi_test(RveNonLinearElastic RVEParameters)
    add_mpi_test(RveNonLinearElasticV2 RVEParameters)
    set(MixedOxideFuelsParameters
        "--mesh"
        "${CMAKE_CURRENT_SOURCE_DIR}/inclusion.msh"
        "--order"
        "1"
        "--library"
        "$<TARGET_FILE:BehaviourTest>"
        "--post-processing"
        "1")
    add_mpi_test(MixedOxideFuels MixedOxideFuelsParameters)
  endif()

endif(MFEM_MGIS_HAVE_TFEL)

# Configure and export example 1 This example is based on UniaxialTensileTest
# from the unit tests

set(EX1_DIR "share/mfem-mgis/examples/ex1")
install(FILES Plasticity.mfront cube.mesh DESTINATION ${EX1_DIR})

install(
  FILES ex1CMakeLists.txt
  DESTINATION ${EX1_DIR}
  RENAME CMakeLists.txt)

install(
  FILES ex1Makefile
  DESTINATION ${EX1_DIR}
  RENAME Makefile)

set(UTILITY_FILE
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/UnitTestingUtilities.hxx
)
add_custom_command(
  OUTPUT ${UTILITY_FILE}
  COMMAND ${PROJECT_SOURCE_DIR}/tests/sed_utility.sh
          ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx ${UTILITY_FILE}
  DEPENDS ${PROJECT_SOURCE_DIR}/tests/sed_utility.sh
          ${PROJECT_SOURCE_DIR}/tests/UnitTestingUtilities.hxx
  COMMENT "seding the UnitTestingUtilities.hxx")

add_custom_target(generate-ex1 ALL DEPENDS ${UTILITY_FILE})

install(FILES ${UTILITY_FILE}
              ${PROJECT_SOURCE_DIR}/tests/UniaxialTensileTest.cxx
        DESTINATION ${EX1_DIR})

set(EX2_DIR "share/mfem-mgis/examples/ex2")
install(FILES PeriodicTest.cxx Elasticity.mfront cube_2mat_per.mesh
        DESTINATION ${EX2_DIR})

install(
  FILES ex2CMakeLists.txt
  DESTINATION ${EX2_DIR}
  RENAME CMakeLists.txt)
