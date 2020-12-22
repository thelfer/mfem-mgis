mfem_mgis_library(MFEMMGIS
  FiniteElementDiscretization.cxx
  PartialQuadratureSpace.cxx
  PartialQuadratureFunction.cxx
  Behaviour.cxx
  Material.cxx
  BehaviourIntegrator.cxx
  BehaviourIntegratorBase.cxx
  BehaviourIntegratorFactory.cxx
  MultiMaterialNonLinearIntegrator.cxx
  ResidualOperator.cxx
  NonLinearEvolutionProblem.cxx
  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.cxx
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.cxx
)

if(enable-openmp)
  target_link_libraries(MFEMMGIS
    PUBLIC OpenMP::OpenMP_CXX)
endif(enable-openmp)