/*!
 * \file   tests/ImposedPressureTest2.cxx
 * \brief  This test tests if the prediction of the solution gives the exact
 * solution for a linear elastic material submited to an external pressure
 * \author Thomas Helfer \date   02/03/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include <iterator>
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
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
       {"NumberOfUniformRefinements", 0},  // faster for testing
       //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", bool(parameters.parallel)}});
  // materials
  problem.addBehaviourIntegrator("Mechanics", 1, parameters.library,
                                 parameters.behaviour);
  auto& m1 = problem.getMaterial(1);
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  const mfem_mgis::real l = 100e9;
  const mfem_mgis::real mu = 75e9;
  mgis::behaviour::setMaterialProperty(m1.s0, "FirstLameCoefficient", l);
  mgis::behaviour::setMaterialProperty(m1.s0, "ShearModulus", mu);
  mgis::behaviour::setMaterialProperty(m1.s1, "FirstLameCoefficient", l);
  mgis::behaviour::setMaterialProperty(m1.s1, "ShearModulus", mu);
  // boundary conditions
  // Only the index is used in this C++ code for manipulating related dof.
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 1, 1));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 2, 2));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 5, 0));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3,
          [](const mfem_mgis::real t) noexcept { return 150e6 * t; }));
  // set the solver parameters
  mfem_mgis::unit_tests::setLinearSolver(problem, parameters);
  problem.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
  // vtk export
  problem.addPostProcessing(
      "ParaviewExportResults",
      {{"OutputFileName",
        "ImposedPressureTest2Output-" + std::string(parameters.behaviour)}});
  // solving the problem in 1 time step
  auto r = problem.solve(0, 1);
  if (!r) {
    return EXIT_FAILURE;
  }
  if (r.iterations != 0) {
    std::cerr << "The newton solver shall not make any iteration";
    return EXIT_FAILURE;
  }
  problem.update();
  mfem_mgis::Profiler::timers::print_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
