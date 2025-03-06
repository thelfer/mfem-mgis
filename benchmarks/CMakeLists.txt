if(MFEM_MGIS_HAVE_TFEL)
	mfront_behaviours_check_library(BehaviourBenchmarks	Elasticity)

	file(COPY
			${PROJECT_SOURCE_DIR}/benchmarks/Elasticity.mfront
			${PROJECT_SOURCE_DIR}/benchmarks/beam-tet.mesh
			DESTINATION ${PROJECT_BINARY_DIR})

# MFEM-MGIS example
	add_executable(MFEMMGISLinearElasticityBenchmark
			MFEMMGISLinearElasticityBenchmark.cxx)
	add_custom_command(
			OUTPUT "src/libBehaviour.so"
			COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Elasticity.mfront
			DEPENDS ${PROJECT_BINARY_DIR}/Elasticity.mfront
			COMMENT "Performing mfront code generation")

  target_link_libraries(MFEMMGISLinearElasticityBenchmark
    PUBLIC MFEMMGIS)


	add_custom_target(generate-libBehaviour-benchmark ALL
			DEPENDS "src/libBehaviour.so")

  add_dependencies(MFEMMGISLinearElasticityBenchmark generate-libBehaviour-benchmark)

# MFEM example
	add_executable(MFEMLinearElasticityBenchmark
			MFEMLinearElasticityBenchmark.cxx)
#  target_include_directories(MFEMLinearElasticityBenchmark
#    PUBLIC ${MFEM_INCLUDE_DIRS})

  target_link_libraries(MFEMLinearElasticityBenchmark
    PUBLIC MFEMMGIS)

endif()
