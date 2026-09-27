/*!
 * \file   tests/tests/MicromorphicDamage2DTest2.cxx
 * \brief
 * \author Thomas Helfer
 * \date   07/12/2021
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "UnitTestingUtilities.hxx"

static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem>
buildMechanicalProblem(
    mgis::Context& ctx,
    const mfem_mgis::unit_tests::TestParameters& test_parameters,
    const mfem_mgis::Parameters& common_problem_parameters) {
  auto or_die = ctx.getFatalFailureHandler();
  constexpr auto E = mfem_mgis::real{200};
  constexpr auto nu = mfem_mgis::real{0.};
  constexpr auto umax = mfem_mgis::real{0.2};
  auto lparameters = common_problem_parameters;
  lparameters.insert(mfem_mgis::throwing, {{"UnknownsSize", 2}});
  auto problem =
      make_shared<mfem_mgis::NonLinearEvolutionProblem>(ctx, lparameters) |
      or_die;
  problem->addBehaviourIntegrator(ctx, "Mechanics", "beam",
                                  test_parameters.library,
                                  "MicromorphicDamageI_SpectralSplit") |
      or_die;
  auto& m = problem->getMaterial(ctx, "beam", 0) | or_die;
  // material properties
  for (const auto& mp : std::map<std::string, double>{{"YoungModulus", E},
                                                      {"PoissonRatio", nu}}) {
    mgis::behaviour::setMaterialProperty(ctx, m.s0, mp.first, mp.second) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s1, mp.first, mp.second) |
        or_die;
  }
  // defining the external state variables
  for (const auto& ev :
       std::map<std::string, double>{{"Temperature", 293.15}, {"Damage", 0}}) {
    mgis::behaviour::setExternalStateVariable(ctx, m.s0, ev.first, ev.second) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s1, ev.first, ev.second) |
        or_die;
  }
  // boundary conditions
  problem->addBoundaryCondition(
      ctx,
      make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          ctx, problem->getFiniteElementDiscretizationPointer(), "left", 0)) |
      or_die;
  problem->addBoundaryCondition(
      ctx,
      make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          ctx, problem->getFiniteElementDiscretizationPointer(), "upper", 1)) |
      or_die;
  problem->addBoundaryCondition(
      ctx,
      make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          ctx, problem->getFiniteElementDiscretizationPointer(), "lower", 1)) |
      or_die;
  problem->addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem->getFiniteElementDiscretizationPointer(), "right",
               0, [](const mfem_mgis::real t) noexcept { return umax * t; }) |
               or_die) |
      or_die;
  // linear solver, convergence critera
  mfem_mgis::unit_tests::setLinearSolver(mfem_mgis::may_abort, ctx, *problem,
                                         test_parameters);
  problem->setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                     {"RelativeTolerance", 1e-4},
                                     {"AbsoluteTolerance", 0},
                                     {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // post-processings
  problem->addPostProcessing(
      ctx, "ParaviewExportResults",
      {{"OutputFileName",
        "MicromorphicDamage2D2TestOutput-MicromorphicDamageI_SpectralSplit"}}) |
      or_die;
  problem->addPostProcessing(
      ctx, "ComputeResultantForceOnBoundary",
      {{"Boundary", "right"},
       {"OutputFileName",
        "MicromorphicDamage2D2TestOutput-"
        "MicromorphicDamageI_SpectralSplit-force.txt"}}) |
      or_die;
  problem->addPostProcessing(
      ctx, "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName",
        "MicromorphicDamage2D2TestIntegrationPointOutput"
        "-MicromorphicDamageI_SpectralSplit"},
       {"Materials", "beam"},
       {"Results",
        std::vector<mfem_mgis::Parameter>{{"EnergyReleaseRate", "Stress"}}}}) |
      or_die;
  return problem;
}

static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem>
buildMicromorphicProblem(
    mgis::Context& ctx,
    const mfem_mgis::unit_tests::TestParameters& test_parameters,
    const mfem_mgis::Parameters& common_problem_parameters) {
  auto or_die = ctx.getFatalFailureHandler();
  constexpr auto Gc = mfem_mgis::real{1};
  constexpr auto l0 = mfem_mgis::real{0.1};
  constexpr auto beta = mfem_mgis::real{300};
  auto lparameters = common_problem_parameters;
  lparameters.insert(mfem_mgis::throwing, {{"UnknownsSize", 1}});
  auto problem =
      std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(ctx, lparameters);
  problem->addBehaviourIntegrator(ctx, "MicromorphicDamage", "beam",
                                  test_parameters.library,
                                  test_parameters.behaviour) |
      or_die;
  auto& m = problem->getMaterial(ctx, "beam", 0) | or_die;
  // material properties
  for (const auto& mp :
       std::map<std::string, double>{{"FractureEnergy", Gc},
                                     {"CharacteristicLength", l0},
                                     {"PenalisationFactor", beta}}) {
    mgis::behaviour::setMaterialProperty(ctx, m.s0, mp.first, mp.second) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, m.s1, mp.first, mp.second) |
        or_die;
  }
  // defining the external state variables
  for (const auto& ev : std::map<std::string, double>{
           {"Temperature", 293.15}, {"EnergyReleaseRate", 0}}) {
    mgis::behaviour::setExternalStateVariable(ctx, m.s0, ev.first, ev.second) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, m.s1, ev.first, ev.second) |
        or_die;
  }
  // linear solver, convergence critera
  mfem_mgis::unit_tests::setLinearSolver(mfem_mgis::may_abort, ctx, *problem,
                                         test_parameters);
  problem->setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                     {"RelativeTolerance", 1e-6},
                                     {"AbsoluteTolerance", 0},
                                     {"MaximumNumberOfIterations", 50}}) |
      or_die;
  // boundary conditions
  problem->addBoundaryCondition(
      ctx,
      make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          ctx, problem->getFiniteElementDiscretizationPointer(), "left", 0) |
          or_die) |
      or_die;
  problem->addBoundaryCondition(
      ctx,
      make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          ctx, problem->getFiniteElementDiscretizationPointer(), "right", 0) |
          or_die) |
      or_die;
  // post-processings
  problem->addPostProcessing(
      ctx, "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "MicromorphicDamage2D2TestIntegrationPointOutput-" +
                              std::string(test_parameters.behaviour)},
       {"Materials", "beam"},
       {"Results",
        std::vector<mfem_mgis::Parameter>{"Damage",
                                          "EnergyReleaseRateValue"}}}) |
      or_die;
  return problem;
}

int main(int argc, char** argv) {
  auto ctx = mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  constexpr auto iter_max = mfem_mgis::size_type{200};
  auto test_parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(test_parameters, argc, argv);
  if (test_parameters.isv_name != nullptr) {
    mfem_mgis::abort("no internal state variable expected");
  }
  //
  const auto common_problem_parameters = mfem_mgis::Parameters{
      {"MeshFileName", test_parameters.mesh_file},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", test_parameters.order},
      {"Hypothesis", "PlaneStrain"},
      {"NumberOfUniformRefinements", test_parameters.parallel ? 1 : 0},
      {"Materials", mfem_mgis::Parameters{{"beam", 5}}},
      {"Boundaries",
       mfem_mgis::Parameters{
           {"left", 3}, {"right", 1}, {"upper", 6}, {"lower", 7}}},
      {"Parallel", bool(test_parameters.parallel)}};
  auto mechanical_problem =
      buildMechanicalProblem(ctx, test_parameters, common_problem_parameters);
  auto micromorphic_problem =
      buildMicromorphicProblem(ctx, test_parameters, common_problem_parameters);
  // solving the problem in 5 time steps, put t1 to 1 and nstep to 100 for the
  // full simulation
  const auto t0 = mfem_mgis::real{0};
  const auto t1 = mfem_mgis::real{0.05};
  const auto nsteps = mfem_mgis::size_type{5};
  const auto dt = (t1 - t0) / nsteps;
  auto t = mfem_mgis::real{0};
  // quadrature functions used to transfer information from one problem to the
  // other
  auto& m1 = micromorphic_problem->getMaterial(ctx, "beam", 0) | or_die;
  mfem_mgis::PartialQuadratureFunction Y(m1.getPartialQuadratureSpacePointer(),
                                         1u);
  mfem_mgis::PartialQuadratureFunction d(m1.getPartialQuadratureSpacePointer(),
                                         1u);
  // using external storage allows to directly modify the values of the
  // quadrature functions Y and d
  auto& m2 = mechanical_problem->getMaterial(ctx, "beam", 0) | or_die;
  mgis::behaviour::setExternalStateVariable(
      ctx, m2.s1, "Damage", d.getValues(),
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) |
      or_die;
  mgis::behaviour::setExternalStateVariable(
      ctx, m1.s1, "EnergyReleaseRate", Y.getValues(),
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) |
      or_die;
  // resolution
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    auto converged = false;
    auto iter = mfem_mgis::size_type{};
    auto mechanical_problem_initial_residual = mfem_mgis::real{};
    auto micromorphic_problem_initial_residual = mfem_mgis::real{};
    std::cout << "\ntime step " << i  //
              << " from " << t << " to " << t + dt << "\n";
    // alternate miminisation algorithm
    while (!converged) {
      std::cout << "time step " << i  //
                << ", alternate minimisation iteration, " << iter << '\n';
      if (iter == 0) {
        mechanical_problem->setSolverParameters(
            ctx, {{"AbsoluteTolerance", 1e-10}}) |
            or_die;
        micromorphic_problem->setSolverParameters(
            ctx, {{"AbsoluteTolerance", 1e-10}}) |
            or_die;
      } else {
        mechanical_problem->setSolverParameters(
            ctx, {{"AbsoluteTolerance",
                   mechanical_problem_initial_residual * 1e-6}}) |
            or_die;
        micromorphic_problem->setSolverParameters(
            ctx, {{"AbsoluteTolerance",
                   micromorphic_problem_initial_residual * 1e-6}}) |
            or_die;
      }
      // solving the mechanical problem
      auto mechanical_output = mechanical_problem->solve(ctx, t, dt);
      if (!mechanical_output.status) {
        mfem_mgis::abort("non convergence of the mechanical problem");
      }
      // passing the energy release rate to the micromorphic problem
      ::mfem_mgis::assign_values(
          ctx, Y,
          mfem_mgis::getInternalStateVariable(ctx, m2, "EnergyReleaseRate") |
              or_die) |
          or_die;
      // solving the micromorphic problem
      auto micromorphic_output = micromorphic_problem->solve(ctx, t, dt);
      if (!micromorphic_output.status) {
        mfem_mgis::abort("non convergence of the micromorphic problem");
      }
      // passing the damage to the mechanical problem
      ::mfem_mgis::assign_values(
          ctx, d,
          mfem_mgis::getInternalStateVariable(ctx, m1, "Damage") | or_die) |
          or_die;
      if (iter == 0) {
        mechanical_problem_initial_residual =
            mechanical_output.initial_residual_norm;
        micromorphic_problem_initial_residual =
            micromorphic_output.initial_residual_norm;
      } else {
        converged = (mechanical_output.iterations == 0) &&
                    (micromorphic_output.iterations == 0);
      }
      ++iter;
      // check convergence
      if ((iter == iter_max) && (!converged)) {
        mfem_mgis::abort("non convergence of the fixed-point problem");
      }
    }
    mechanical_problem->executePostProcessings(ctx, t, dt) | or_die;
    micromorphic_problem->executePostProcessings(ctx, t, dt) | or_die;
    mechanical_problem->update(ctx) | or_die;
    micromorphic_problem->update(ctx) | or_die;
    t += dt;
  }
  return EXIT_SUCCESS;
}
