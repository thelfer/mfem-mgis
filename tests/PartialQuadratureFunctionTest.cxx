/*!
 * \file   tests/UniaxialTensileTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
#ifdef DO_USE_MPI
  static constexpr const auto parallel = true;
#else
  static constexpr const auto parallel = false;
#endif
  constexpr const auto dim = mfem_mgis::size_type{3};
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  auto success = true;
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", parameters.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", parameters.order},
       {"UnknownsSize", dim},
       {"NumberOfUniformRefinements", parallel ? 2 : 0},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", parallel}});
  // materials
  problem.addBehaviourIntegrator("Mechanics", 1, parameters.library,
                                 parameters.behaviour);
  auto& m1 = problem.getMaterial(1);
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  //
  auto e = getGradient(m1, "Strain");
  auto s = getThermodynamicForce(m1, "Stress");
  auto p = getInternalStateVariable(m1, "EquivalentStrain");
  success = e.getNumberOfComponents() >= 0 && s.getNumberOfComponents() >= 0 &&
            p.getNumberOfComponents() >= 0;
  ;

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
