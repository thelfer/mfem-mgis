/*!
 * \file   BehaviourIntegratorsSelectionTest.cxx
 * \brief  Test of the selection of the behaviour integrators used by the
 * post-processings
 * \date   02/10/2026
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto params = unit_tests::TestParameters{};
  unit_tests::parseCommandLineOptions(params, argc, argv);
  auto ctx = Context{};
  auto or_die = ctx.getFatalFailureHandler();
  auto problem = construct<NonLinearEvolutionProblem>(
                     ctx, Parameters{{"MeshFileName", params.mesh_file},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", bool(params.parallel)}}) |
                 or_die;
  // two identical behaviour integrators on each material
  for (const size_type m : {1, 2}) {
    for (size_type b = 0; b != 2; ++b) {
      problem.addBehaviourIntegrator(ctx, "Mechanics", m, params.library,
                                     params.behaviour) |
          or_die;
      auto& mat = problem.getMaterial(ctx, m, b) | or_die;
      for (auto* s : {&mat.s0, &mat.s1}) {
        mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature",
                                                  293.15) |
            or_die;
        mgis::behaviour::setMaterialProperty(ctx, *s, "YoungModulus", 200e9) |
            or_die;
        mgis::behaviour::setMaterialProperty(ctx, *s, "PoissonRatio", 0.3) |
            or_die;
      }
    }
  }
  // clamped on the first boundary, traction along x on the second one
  auto pfed = problem.getFiniteElementDiscretizationPointer();
  for (size_type c = 0; c != 3; ++c) {
    problem.addBoundaryCondition(
        ctx, make_unique<UniformDirichletBoundaryCondition>(ctx, pfed, 1, c) |
                 or_die) |
        or_die;
  }
  problem.addBoundaryCondition(
      ctx,
      make_unique<UniformDirichletBoundaryCondition>(
          ctx, pfed, 2, 0, [](const real t) noexcept { return 1e-3 * t; }) |
          or_die) |
      or_die;
  unit_tests::setLinearSolver(may_abort, ctx, problem, params);
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // post-processings using all the behaviour integrators
  for (const auto* const n : {"ComputeResultantForceOnBoundary",
                              "MeanThermodynamicForces", "StoredEnergy"}) {
    auto pp = Parameters{
        {"OutputFileName", std::string{n} + "-BehaviourIntegratorsSelection-" +
                               std::to_string(params.parallel) + ".txt"},
        {"BehaviourIntegrator", "All"}};
    if (std::string_view{n} == "ComputeResultantForceOnBoundary") {
      pp.insert(throwing, {{"Boundary", 1}});
    }
    problem.addPostProcessing(ctx, n, pp) | or_die;
  }
  if (!problem.solve(ctx, 0, 1).status) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  if (!problem.executePostProcessings(ctx, 0, 1)) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return EXIT_FAILURE;
  }
  // first component of the contribution of the process to the resultant on
  // the first boundary
  const auto edofs = getElementsDegreesOfFreedomOnBoundary(problem, 1);
  auto resultant =
      [&](const BehaviourIntegratorsSelection& s) -> std::optional<real> {
    auto lctx = Context{};
    auto F = mfem::Vector{};
    if (!computeResultantForceOnBoundary(lctx, F, problem, edofs, s)) {
      return {};
    }
    return F[0];
  };
  const auto F0 = resultant({.index = 0});
  const auto F1 = resultant({.index = 1});
  const auto Fa = resultant({.all = true});
  auto success = true;
  // the selection is ambiguous without the `BehaviourIntegrator` parameter
  success = !resultant({}).has_value() && success;
  // no third behaviour integrator
  success = !resultant({.index = 2}).has_value() && success;
  success = F0.has_value() && F1.has_value() && Fa.has_value() && success;
  if (success) {
    const auto eps = 1e-12 * std::abs(*Fa);
    success = (std::abs(*F0 - *F1) <= eps) &&
              (std::abs(*Fa - (*F0 + *F1)) <= eps) && success;
  }
  if (!success) {
    std::cerr << "invalid selection of the behaviour integrators\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
