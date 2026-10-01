/*!
 * \file   tests/MicromorphicDamage2DTest.cxx
 * \brief  Tests of the bidimensional micromorphic damage behaviour integrator
 * \author Thomas Helfer
 * \date   07/12/2021
 *
 * This test compares the solution obtained with the MicromorphicDamage
 * behaviour integrator in 2D on a bar with an analytical solution using
 * the MiehePhaseFieldDamage behaviour. A specified history function is
 * imposed so that the solution of the phase field equation is:
 * \f[
 * d(x,y) = \sin(4\,\pi\,x)/2
 * \f]
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  constexpr auto pi = mfem_mgis::real{3.14159265358979323846};
  constexpr auto Gc = mfem_mgis::real{1};
  constexpr auto l0 = mfem_mgis::real{0.1};
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  if (parameters.isv_name != nullptr) {
    mfem_mgis::abort("no internal state variable expected");
  }
  // building the non linear problem
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx, mfem_mgis::Parameters{{"MeshFileName", parameters.mesh_file},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", parameters.order},
                                     {"UnknownsSize", 1},
                                     {"NumberOfUniformRefinements",
                                      parameters.parallel ? 1 : 0},
                                     {"Hypothesis", "PlaneStrain"},
                                     {"Parallel", bool(parameters.parallel)}}) |
      or_die;
  // materials
  problem.addBehaviourIntegrator(ctx, "MicromorphicDamage", 5,
                                 parameters.library, parameters.behaviour) |
      or_die;
  auto& m = problem.getMaterial(ctx, 5, 0) | or_die;
  for (const auto& ev : std::map<std::string, double>{{"Temperature", 293.15},
                                                      {"HistoryFunction", 0}}) {
    mgis::behaviour::setExternalStateVariable(ctx, m.s0, ev.first, ev.second) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s1, ev.first, ev.second) |
        or_die;
  }
  //
  const auto H = mfem_mgis::PartialQuadratureFunction::evaluate(
      m.getPartialQuadratureSpacePointer(),
      [](const mfem_mgis::real x, const mfem_mgis::real) noexcept {
        const auto d = std::sin(4 * pi * x) / 2;
        return Gc * d * (1 + l0 * l0 * 16 * pi * pi) / (2 * l0 * (1 - d));
      });
  mgis::behaviour::setExternalStateVariable(
      ctx, m.s1, "HistoryFunction", H->getValues(),
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) |
      or_die;
  // material properties
  for (const auto& mp : std::map<std::string, double>{
           {"FractureEnergy", Gc}, {"RegularizationLength", l0}}) {
    mgis::behaviour::setMaterialProperty(ctx, m.s0, mp.first, mp.second) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s1, mp.first, mp.second) |
        or_die;
  }
  // boundary conditions
  //    $PhysicalNames
  // 1 3 "LD"
  // 1 1 "LG"
  //    $EndPhysicalNames
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3, 0) |
               or_die) |
      or_die;
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 1, 0) |
               or_die) |
      or_die;
  // set the solver parameters
  mfem_mgis::unit_tests::setLinearSolver(mfem_mgis::may_abort, ctx, problem,
                                         parameters);
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-10},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // vtk export
  problem.addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName", "MicromorphicDamage2DTestOutput-" +
                              std::string(parameters.behaviour)}}) |
      or_die;
  // solving the problem in 1 time steps
  const auto t0 = mfem_mgis::real{0};
  const auto t1 = mfem_mgis::real{1};
  const auto nsteps = mfem_mgis::size_type{1};
  const auto dt = (t1 - t0) / nsteps;
  auto t = mfem_mgis::real{0};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    problem.solve(ctx, t, dt) | or_die;
    problem.executePostProcessings(ctx, t, dt) | or_die;
    problem.update(ctx) | or_die;
    t += dt;
  }
  //
  const auto osuccess = mfem_mgis::compareToAnalyticalSolution(
      ctx, problem,
      [](mfem::Vector& u, const mfem::Vector& x) {
        u[0] = std::sin(4 * pi * x[0]) / 2;
      },
      {{"CriterionThreshold", 3e-6}});
  if (mfem_mgis::isInvalid(osuccess)) {
    mfem_mgis::getErrorStream() << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  return *osuccess ? EXIT_SUCCESS : EXIT_FAILURE;
}
