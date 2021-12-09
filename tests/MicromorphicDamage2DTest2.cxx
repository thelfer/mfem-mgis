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
    const mfem_mgis::unit_tests::TestParameters& test_parameters,
    const mfem_mgis::Parameters& common_problem_parameters) {
  constexpr auto E = mfem_mgis::real{200};
  constexpr auto nu = mfem_mgis::real{0.};
  constexpr auto umax = mfem_mgis::real{0.2};
  auto lparameters = common_problem_parameters;
  lparameters.insert({{"UnknownsSize", 2}});
  auto problem =
      std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(lparameters);
  problem->addBehaviourIntegrator("Mechanics", "beam", test_parameters.library,
                                  "MicromorphicDamageI_SpectralSplit");
  auto& m = problem->getMaterial("beam");
  // material properties
  for (const auto& mp : std::map<std::string, double>{{"YoungModulus", E},
                                                      {"PoissonRatio", nu}}) {
    mgis::behaviour::setMaterialProperty(m.s0, mp.first, mp.second);
    mgis::behaviour::setMaterialProperty(m.s1, mp.first, mp.second);
  }
  // defining the external state variables
  for (const auto& ev :
       std::map<std::string, double>{{"Temperature", 293.15}, {"Damage", 0}}) {
    mgis::behaviour::setExternalStateVariable(m.s0, ev.first, ev.second);
    mgis::behaviour::setExternalStateVariable(m.s1, ev.first, ev.second);
  }
  // boundary conditions
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "left", 0));
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "upper",
          1));
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "lower",
          1));
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "right", 0,
          [](const mfem_mgis::real t) { return umax * t; }));
  // linear solver, convergence critera
  mfem_mgis::unit_tests::setLinearSolver(*problem, test_parameters);
  problem->setSolverParameters({{"VerbosityLevel", 0},
                                {"RelativeTolerance", 1e-4},
                                {"AbsoluteTolerance", 0},
                                {"MaximumNumberOfIterations", 10}});
  // post-processings
  problem->addPostProcessing(
      "ParaviewExportResults",
      {{"OutputFileName",
        "MicromorphicDamage2D2TestOutput-MicromorphicDamageI_SpectralSplit"}});
  problem->addPostProcessing("ComputeResultantForceOnBoundary",
                             {{"Boundary", "right"},
                              {"OutputFileName",
                               "MicromorphicDamage2D2TestOutput-"
                               "MicromorphicDamageI_SpectralSplit-force.txt"}});
  problem->addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName",
        "MicromorphicDamage2D2TestIntegrationPointOutput"
        "-MicromorphicDamageI_SpectralSplit"},
       {"Materials", "beam"},
       {"Results",
        std::vector<mfem_mgis::Parameter>{{"EnergyReleaseRate", "Stress"}}}});
  return problem;
}

static std::shared_ptr<mfem_mgis::NonLinearEvolutionProblem>
buildMicromorphicProblem(
    const mfem_mgis::unit_tests::TestParameters& test_parameters,
    const mfem_mgis::Parameters& common_problem_parameters) {
  constexpr auto Gc = mfem_mgis::real{1};
  constexpr auto l0 = mfem_mgis::real{0.1};
  constexpr auto beta = mfem_mgis::real{300};
  auto lparameters = common_problem_parameters;
  lparameters.insert({{"UnknownsSize", 1}});
  auto problem =
      std::make_shared<mfem_mgis::NonLinearEvolutionProblem>(lparameters);
  problem->addBehaviourIntegrator("MicromorphicDamage", "beam",
                                  test_parameters.library,
                                  test_parameters.behaviour);
  auto& m = problem->getMaterial("beam");
  // material properties
  for (const auto& mp :
       std::map<std::string, double>{{"FractureEnergy", Gc},
                                     {"CharacteristicLength", l0},
                                     {"PenalisationFactor", beta}}) {
    mgis::behaviour::setMaterialProperty(m.s0, mp.first, mp.second);
    mgis::behaviour::setMaterialProperty(m.s1, mp.first, mp.second);
  }
  // defining the external state variables
  for (const auto& ev : std::map<std::string, double>{
           {"Temperature", 293.15}, {"EnergyReleaseRate", 0}}) {
    mgis::behaviour::setExternalStateVariable(m.s0, ev.first, ev.second);
    mgis::behaviour::setExternalStateVariable(m.s1, ev.first, ev.second);
  }
  // linear solver, convergence critera
  mfem_mgis::unit_tests::setLinearSolver(*problem, test_parameters);
  problem->setSolverParameters({{"VerbosityLevel", 0},
                                {"RelativeTolerance", 1e-6},
                                {"AbsoluteTolerance", 0},
                                {"MaximumNumberOfIterations", 50}});
  // boundary conditions
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "left", 0));
  problem->addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem->getFiniteElementDiscretizationPointer(), "right", 0));
  // post-processings
  problem->addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "MicromorphicDamage2D2TestIntegrationPointOutput-" +
                              std::string(test_parameters.behaviour)},
       {"Materials", "beam"},
       {"Results",
        std::vector<mfem_mgis::Parameter>{"Damage", "EnergyReleaseRateValue"}}});
  return problem;
}

static void extractInternalStateVariable(
    mfem_mgis::PartialQuadratureFunction& f,
    mgis::behaviour::MaterialStateManager& s,
    const mgis::string_view n) {
  // checking compatibility
  const auto& q = f.getPartialQuadratureSpace();
  if (q.getNumberOfIntegrationPoints() !=
      static_cast<mfem_mgis::size_type>(s.n)) {
    mfem_mgis::raise(
        "extractInternalStateVariable: "
        "unmatched number of integration points");
  }
  const auto& iv = mgis::behaviour::getVariable(s.b.isvs, n);
  const auto nc = mgis::behaviour::getVariableSize(iv, s.b.hypothesis);
  if (f.getNumberOfComponents() != nc) {
    mfem_mgis::raise(
        "extractInternalStateVariable: "
        "unmatched number of components");
  }
  //
  if (nc != 1u) {
    mfem_mgis::raise(
        "extractInternalStateVariable: "
        "unsupported number of components");
  }
  //
  const auto o =
      mgis::behaviour::getVariableOffset(s.b.isvs, n, s.b.hypothesis);
  const auto stride = s.internal_state_variables_stride;
  auto* const p = f.getValues().data();
  const auto* const piv = s.internal_state_variables.data() + o;
  for (mfem_mgis::size_type i = 0; i != q.getNumberOfIntegrationPoints(); ++i) {
    p[i] = piv[i * stride];
  }
}  // end of extractInternalStateVariable

int main(int argc, char** argv) {
  constexpr auto iter_max = mfem_mgis::size_type{200};
#ifdef DO_USE_MPI
  static constexpr const auto parallel = true;
#else
  static constexpr const auto parallel = false;
#endif
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
      {"NumberOfUniformRefinements", parallel ? 2 : 0},
      {"Materials", mfem_mgis::Parameters{{"beam", 5}}},
      {"Boundaries",
       mfem_mgis::Parameters{
           {"left", 3}, {"right", 1}, {"upper", 6}, {"lower", 7}}},
      {"Parallel", parallel}};
  auto mechanical_problem =
      buildMechanicalProblem(test_parameters, common_problem_parameters);
  auto micromorphic_problem =
      buildMicromorphicProblem(test_parameters, common_problem_parameters);
  // solving the problem in 100 time steps
  const auto t0 = mfem_mgis::real{0};
  const auto t1 = mfem_mgis::real{1};
  const auto nsteps = mfem_mgis::size_type{100};
  const auto dt = (t1 - t0) / nsteps;
  auto t = mfem_mgis::real{0};
  // quadrature functions used to transfer information from one problem to the
  // other
  mfem_mgis::PartialQuadratureFunction Y(
      micromorphic_problem->getMaterial("beam")
          .getPartialQuadratureSpacePointer(),
      1u);
  mfem_mgis::PartialQuadratureFunction d(
      micromorphic_problem->getMaterial("beam")
          .getPartialQuadratureSpacePointer(),
      1u);
  // using external storage allows to directly modify the values of the
  // quadrature functions Y and d
  mgis::behaviour::setExternalStateVariable(
      mechanical_problem->getMaterial("beam").s1, "Damage", d.getValues(),
      mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE);
  mgis::behaviour::setExternalStateVariable(
      micromorphic_problem->getMaterial("beam").s1, "EnergyReleaseRate",
      Y.getValues(), mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE);
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
        mechanical_problem->setSolverParameters({{"AbsoluteTolerance", 1e-10}});
        micromorphic_problem->setSolverParameters(
            {{"AbsoluteTolerance", 1e-10}});
      } else {
        mechanical_problem->setSolverParameters(
            {{"AbsoluteTolerance",
              mechanical_problem_initial_residual * 1e-6}});
        micromorphic_problem->setSolverParameters(
            {{"AbsoluteTolerance",
              micromorphic_problem_initial_residual * 1e-6}});
      }
      // solving the mechanical problem
      auto mechanical_output = mechanical_problem->solve(t, dt);
      if (!mechanical_output.status) {
        mfem_mgis::raise("non convergence of the mechanical problem");
      }
      // passing the energy release rate to the micromorphic problem
      extractInternalStateVariable(
          Y, mechanical_problem->getMaterial("beam").s1, "EnergyReleaseRate");
      // solving the micromorphic problem
      auto micromorphic_output = micromorphic_problem->solve(t, dt);
      if (!micromorphic_output.status) {
        mfem_mgis::raise("non convergence of the micromorphic problem");
      }
      // passing the damage to the mechanical problem
      extractInternalStateVariable(
          d, micromorphic_problem->getMaterial("beam").s1, "Damage");
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
        mfem_mgis::raise("non convergence of the fixed-point problem");
      }
    }
    mechanical_problem->executePostProcessings(t, dt);
    micromorphic_problem->executePostProcessings(t, dt);
    mechanical_problem->update();
    micromorphic_problem->update();
    t += dt;
  }
  return EXIT_SUCCESS;
}
