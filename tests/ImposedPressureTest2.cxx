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
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  constexpr const auto dim = mfem_mgis::size_type{3};
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  auto success = true;
  // building the non linear problem
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", parameters.mesh_file},
              {"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", parameters.order},
              {"UnknownsSize", dim},
              {"NumberOfUniformRefinements", 0},  // faster for testing
              //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
              {"Hypothesis", "Tridimensional"},
              {"Parallel", bool(parameters.parallel)}}) |
      or_die;
  // materials
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, parameters.library,
                                 parameters.behaviour) |
      or_die;
  auto& m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s0, "Temperature", 293.15) |
      or_die;
  mgis::behaviour::setExternalStateVariable(ctx, m1.s1, "Temperature", 293.15) |
      or_die;
  const mfem_mgis::real l = 100e9;
  const mfem_mgis::real mu = 75e9;
  mgis::behaviour::setMaterialProperty(ctx, m1.s0, "FirstLameCoefficient", l) |
      or_die;
  mgis::behaviour::setMaterialProperty(ctx, m1.s0, "ShearModulus", mu) | or_die;
  mgis::behaviour::setMaterialProperty(ctx, m1.s1, "FirstLameCoefficient", l) |
      or_die;
  mgis::behaviour::setMaterialProperty(ctx, m1.s1, "ShearModulus", mu) | or_die;
  // boundary conditions
  // Only the index is used in this C++ code for manipulating related dof.
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 1, 1) |
               or_die) |
      or_die;
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 2, 2) |
               or_die) |
      or_die;
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 5, 0) |
               or_die) |
      or_die;
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3,
               [](const mfem_mgis::real t) noexcept { return 150e6 * t; }) |
               or_die) |
      or_die;
  // set the solver parameters
  mfem_mgis::unit_tests::setLinearSolver(mfem_mgis::may_abort, ctx, problem,
                                         parameters);
  problem.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // vtk export
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName",
        "ImposedPressureTest2Output-" + std::string(parameters.behaviour)}}) |
      or_die;
  // solving the problem in 1 time step
  auto r = problem.solve(ctx, 0, 1);
  if (!r) {
    return EXIT_FAILURE;
  }
  if (r.iterations != 0) {
    std::cerr << "The newton solver shall not make any iteration";
    return EXIT_FAILURE;
  }
  problem.update(ctx) | or_die;
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
