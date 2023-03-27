/*!
 * \file   tests/MicromorphicDamage2DTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   07/12/2021
 *
 * This test compares the solution obtained with the MicromorphicDamage
 * behaviour integrator in 2D on a bar with an analytical solution using
 * the MiehePhaseFieldDamage behaviour. A specified history function is
 * imposed so that the solution of the phase field equation is:
 * \f[
 * d(x,y) = \sin(4\,\,\pi\,x)/2
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
  constexpr auto pi = mfem_mgis::real{3.14159265358979323846};
  constexpr auto Gc = mfem_mgis::real{1};
  constexpr auto l0 = mfem_mgis::real{0.1};
#ifdef DO_USE_MPI
  static constexpr const auto parallel = true;
#else
  static constexpr const auto parallel = false;
#endif
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  if (parameters.isv_name != nullptr) {
    mfem_mgis::abort("no internal state variable expected");
  }
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", parameters.mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", parameters.order},
       {"UnknownsSize", 1},
       {"NumberOfUniformRefinements", parallel ? 2 : 0},
       {"Hypothesis", "PlaneStrain"},
       {"Parallel", parallel}});
  // materials
  problem.addBehaviourIntegrator("MicromorphicDamage", 5, parameters.library,
                                 parameters.behaviour);
  auto& m = problem.getMaterial(5);
  for (const auto& ev : std::map<std::string, double>{{"Temperature", 293.15},
                                                      {"HistoryFunction", 0}}) {
    mgis::behaviour::setExternalStateVariable(m.s0, ev.first, ev.second);
    mgis::behaviour::setExternalStateVariable(m.s1, ev.first, ev.second);
  }
  //
  const auto H = mfem_mgis::PartialQuadratureFunction::evaluate(
      m.getPartialQuadratureSpacePointer(),
      [](const mfem_mgis::real x, const mfem_mgis::real) {
        const auto d = std::sin(4 * pi * x) / 2;
        return Gc * d * (1 + l0 * l0 * 16 * pi * pi) / (2 * l0 * (1 - d));
      });
  mgis::behaviour::setExternalStateVariable(
      m.s1, "HistoryFunction", H->getValues(),
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE);
  // material properties
  for (const auto& mp : std::map<std::string, double>{
           {"FractureEnergy", Gc}, {"RegularizationLength", l0}}) {
    mgis::behaviour::setMaterialProperty(m.s0, mp.first, mp.second);
    mgis::behaviour::setMaterialProperty(m.s1, mp.first, mp.second);
  }
  // boundary conditions
  //    $PhysicalNames
  // 1 3 "LD"
  // 1 1 "LG"
  //    $EndPhysicalNames
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3, 0));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 1, 0));
  // set the solver parameters
  mfem_mgis::unit_tests::setLinearSolver(problem, parameters);
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
  // vtk export
  problem.addPostProcessing(
      "ParaviewExportResults",
      {{"OutputFileName", "MicromorphicDamage2DTestOutput-" +
                              std::string(parameters.behaviour)}});
  // solving the problem in 1 time steps
  const auto t0 = mfem_mgis::real{0};
  const auto t1 = mfem_mgis::real{1};
  const auto nsteps = mfem_mgis::size_type{1};
  const auto dt = (t1 - t0) / nsteps;
  auto t = mfem_mgis::real{0};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    if (!problem.solve(t, dt)) {
      mfem_mgis::raise("non convergence");
    }
    problem.executePostProcessings(t, dt);
    problem.update();
    t += dt;
  }
  //
  const auto success = mfem_mgis::compareToAnalyticalSolution(
      problem,
      [](mfem::Vector& u, const mfem::Vector& x) {
        u[0] = std::sin(4 * pi * x[0]) / 2;
      },
      {{"CriterionThreshold", 3e-6}});
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
