/*!
 * \file   tests/FiniteElementSpaceOrderingTest.cxx
 * \brief  This test checks that the results do not depend on the ordering of
 * the degrees of freedom (`FiniteElementSpaceOrdering` parameter).
 *
 * A cube, clamped on the face x = 0, is loaded by a pressure on the face
 * x = 1. The problem is solved with the degrees of freedom ordered by nodes
 * and by vector dimension. The thermodynamic forces at the integration points
 * and the resultant force on the clamped face must be the same.
 *
 * \date   05/10/2026
 */

#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int parallel = 0;
};

//! \brief results compared between the two orderings
struct Results {
  //! \brief thermodynamic forces at the integration points
  std::vector<mfem_mgis::real> thermodynamic_forces;
  //! \brief resultant force on the clamped face, only read by the main process
  std::array<mfem_mgis::real, 3> resultant_force = {0, 0, 0};
};

//! \brief imposed pressure
constexpr auto pressure = mfem_mgis::real{150e6};

static TestParameters parseCommandLineOptions(int argc, char** argv) {
  auto p = TestParameters{};
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (p.mesh_file == nullptr) || (p.library == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  return p;
}

/*!
 * \return the resultant force of the last line of the given file
 * \param[in] f: file written by the `ComputeResultantForceOnBoundary`
 * post-processing
 */
static std::array<mfem_mgis::real, 3> readResultantForce(const std::string& f) {
  auto in = std::ifstream(f);
  auto line = std::string{};
  auto last = std::string{};
  while (std::getline(in, line)) {
    if ((!line.empty()) && (line[0] != '#')) {
      last = line;
    }
  }
  auto t = mfem_mgis::real{};
  auto F = std::array<mfem_mgis::real, 3>{};
  if (!(std::istringstream(last) >> t >> F[0] >> F[1] >> F[2])) {
    mfem_mgis::raise("can't read the resultant force in '" + f + "'");
  }
  return F;
}

/*!
 * \return the results obtained with the given ordering
 * \param[in, out] ctx: execution context
 * \param[in] mesh: mesh discretization
 * \param[in] p: test parameters
 * \param[in] ordering: ordering of the degrees of freedom
 */
static Results solve(mfem_mgis::Context& ctx,
                     mfem_mgis::MeshDiscretization& mesh,
                     const TestParameters& p,
                     const std::string& ordering) {
  auto or_die = ctx.getFatalFailureHandler();
  auto problem =
      mfem_mgis::construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx, mesh,
          mfem_mgis::Parameters{{"FiniteElementFamily", "H1"},
                                {"FiniteElementOrder", 2},
                                {"FiniteElementSpaceOrdering", ordering},
                                {"UnknownsSize", 3},
                                {"Hypothesis", "Tridimensional"}}) |
      or_die;
  // material
  problem.addBehaviourIntegrator(ctx, "Mechanics", "Cube", p.library,
                                 "IsotropicLinearElasticity") |
      or_die;
  auto& m = problem.getMaterial(ctx, "Cube", 0) | or_die;
  for (auto* const s : {&(m.s0), &(m.s1)}) {
    mgis::behaviour::setExternalStateVariable(ctx, *s, "Temperature", 293.15) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, *s, "FirstLameCoefficient",
                                         100e9) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, *s, "ShearModulus", 75e9) |
        or_die;
  }
  // boundary conditions
  for (int c = 0; c != 3; ++c) {
    problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", "Xmin"}, {"Component", c}}) |
        or_die;
  }
  problem.addBoundaryCondition(
      ctx, make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), "Xmax",
               [](const mfem_mgis::real t) noexcept { return pressure * t; }) |
               or_die) |
      or_die;
  // solver parameters
  problem.setLinearSolver(ctx, "CGSolver",
                          {{"VerbosityLevel", 0},
                           {"AbsoluteTolerance", 0.},
                           {"RelativeTolerance", 1e-14},
                           {"MaximumNumberOfIterations", 5000}}) |
      or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 0},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  // resultant force on the clamped face
  const auto f = "FiniteElementSpaceOrderingTest-" + ordering +
                 (p.parallel ? "-parallel" : "") + ".txt";
  problem.addPostProcessing(ctx, "ComputeResultantForceOnBoundary",
                            {{"Boundary", "Xmin"}, {"OutputFileName", f}}) |
      or_die;
  // resolution
  problem.solve(ctx, 0, 1) | or_die;
  problem.executePostProcessings(ctx, 0, 1) | or_die;
  auto r = Results{};
  r.thermodynamic_forces.assign(m.s1.thermodynamic_forces.begin(),
                                m.s1.thermodynamic_forces.end());
  // the file is written by the main process
  if (mfem_mgis::isMainProcess(mesh)) {
    r.resultant_force = readResultantForce(f);
  }
  return r;
}

int main(int argc, char** argv) {
  auto ctx = mfem_mgis::Context{};
  mfem_mgis::initialize(argc, argv);
  const auto p = parseCommandLineOptions(argc, argv);
  auto or_die = ctx.getFatalFailureHandler();
  // both problems are built on the same mesh
  auto mesh =
      mfem_mgis::construct<mfem_mgis::MeshDiscretization>(
          ctx,
          mfem_mgis::Parameters{
              {"MeshFileName", p.mesh_file},
              {"NumberOfUniformRefinements", 2},
              {"Materials", mfem_mgis::Parameters{{"Cube", 1}}},
              {"Boundaries", mfem_mgis::Parameters{{"Xmax", 3}, {"Xmin", 5}}},
              {"Parallel", bool(p.parallel)}}) |
      or_die;
  const auto r1 = solve(ctx, mesh, p, "byNODES");
  const auto r2 = solve(ctx, mesh, p, "byVDIM");
  // absolute tolerance on the stresses and on the resultant force (the area
  // of the loaded face is 1)
  constexpr auto eps = 1e-8 * pressure;
  auto success =
      r1.thermodynamic_forces.size() == r2.thermodynamic_forces.size();
  if (success) {
    for (std::size_t i = 0; i != r1.thermodynamic_forces.size(); ++i) {
      const auto e =
          std::abs(r1.thermodynamic_forces[i] - r2.thermodynamic_forces[i]);
      if (e > eps) {
        mfem_mgis::getErrorStream()
            << "invalid thermodynamic force " << i << ": "
            << r1.thermodynamic_forces[i] << " (byNODES) vs "
            << r2.thermodynamic_forces[i] << " (byVDIM)\n";
        success = false;
      }
    }
  } else {
    mfem_mgis::getErrorStream() << "unmatched number of integration points\n";
  }
  // the resultant force balances the pressure in both cases
  if (mfem_mgis::isMainProcess(mesh)) {
    const auto expected = std::array<mfem_mgis::real, 3>{pressure, 0, 0};
    for (const auto& [n, F] : {std::pair{"byNODES", r1.resultant_force},
                               std::pair{"byVDIM", r2.resultant_force}}) {
      for (std::size_t i = 0; i != 3; ++i) {
        if (std::abs(std::abs(F[i]) - expected[i]) > eps) {
          mfem_mgis::getErrorStream()
              << "invalid resultant force (" << n << "), component " << i
              << ": " << F[i] << '\n';
          success = false;
        }
      }
    }
  }
  success = mfem_mgis::isTrueOnAllProcesses(mesh, success);
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
