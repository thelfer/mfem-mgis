mfem_mgis_library(MFEMMGIS
  PartialQuadratureSpace.cxx
  PartialQuadratureFunction.cxx
  Behaviour.cxx
  Material.cxx
  BehaviourIntegrator.cxx
  BehaviourIntegratorBase.cxx
  SmallStrainMechanicalBehaviourIntegrator.cxx
  BehaviourIntegratorFactory.cxx
  MGISIntegrator.cxx
  ResidualOperator.cxx
  NonLinearEvolutionProblem.cxx
)

if(enable-openmp)
  target_link_libraries(MFEMMGIS
    PUBLIC OpenMP::OpenMP_CXX)
endif(enable-openmp)