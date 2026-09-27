#include <cstdlib>
#include "mfem/linalg/sparsemat.hpp"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/lininteg.hpp"
#include "mfem/fem/datacollection.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pbilinearform.hpp"
#endif /* MFEM_USE_MPI */
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/LinearSolverHandler.hxx"
#include "MFEMMGIS/AbstractDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(int argc, char *argv[]) {
  //
  mfem_mgis::initialize(argc, argv);
  //
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  // parse command-line options.
  const char *mesh_file = nullptr;
  const char *library = nullptr;
  int order = 1;
  int parallel = 0;

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library",
                 "library containing the behaviour.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    return EXIT_FAILURE;
  }
  if (mesh_file == nullptr) {
    mfem_mgis::getOutputStream() << "no mesh file specified\n";
    args.PrintUsage(mfem_mgis::getOutputStream());
    return EXIT_FAILURE;
  }
  if (library == nullptr) {
    mfem_mgis::getOutputStream() << "no library specified\n";
    args.PrintUsage(mfem_mgis::getOutputStream());
    return EXIT_FAILURE;
  }
  args.PrintOptions(mfem_mgis::getOutputStream());
  //
  auto problem =
      construct<mfem_mgis::NonLinearEvolutionProblem>(
          ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", order},
                                     {"UnknownsSize", 3},
                                     {"NumberOfUniformRefinements", 2},
                                     {"Hypothesis", "Tridimensional"},
                                     {"Parallel", bool(parallel)}}) |
      or_die;
  //
  problem.addBehaviourIntegrator(ctx, "Mechanics", 1, library, "Elasticity") |
      or_die;
  auto &m1 = problem.getMaterial(ctx, 1, 0) | or_die;
  for (auto *ps : {&m1.s0, &m1.s1}) {
    mgis::behaviour::setMaterialProperty(ctx, *ps, "FirstLameCoefficient",
                                         100e9) |
        or_die;
    mgis::behaviour::setMaterialProperty(ctx, *ps, "ShearModulus", 75e9) |
        or_die;
    mgis::behaviour::setExternalStateVariable(ctx, *ps, "Temperature", 293.15) |
        or_die;
  }
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
      ctx, make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
               ctx, problem.getFiniteElementDiscretizationPointer(), 3, 0,
               [](const auto t) noexcept { return 3e-2 * t; }) |
               or_die) |
      or_die;
  problem.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
  // set the solver parameters
  problem.setLinearSolver(ctx, "CGSolver",
                          {{"VerbosityLevel", 1},
                           {"AbsoluteTolerance", 1e-16},
                           {"RelativeTolerance", 1e-16},
                           {"MaximumNumberOfIterations", 1000}}) |
      or_die;
  problem.setSolverParameters(ctx, {{"VerbosityLevel", 2},
                                    {"RelativeTolerance", 1e-12},
                                    {"AbsoluteTolerance", 0.},
                                    {"MaximumNumberOfIterations", 10}}) |
      or_die;
  //
  auto r = problem.solve(ctx, 0, 1);
  if (!r) {
    std::cout << "Non convergence of the nonlinear algorithm\n";
    return EXIT_FAILURE;
  }
  if (r.iterations != 0) {
    std::cerr << "The newton solver shall not make any iteration";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
