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
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(int argc, char *argv[]) {
  //
  mfem_mgis::initialize(argc, argv);
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
  auto ctx = mfem_mgis::Context{};
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", 3},
       {"NumberOfUniformRefinements", 2},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", bool(parallel)}});
  //
  problem.addBehaviourIntegrator("Mechanics", 1, library, "Elasticity");
  auto &m1 = problem.getMaterial(1);
  for (auto *ps : {&m1.s0, &m1.s1}) {
    mgis::behaviour::setMaterialProperty(*ps, "FirstLameCoefficient", 100e9);
    mgis::behaviour::setMaterialProperty(*ps, "ShearModulus", 75e9);
    mgis::behaviour::setExternalStateVariable(*ps, "Temperature", 293.15);
  }
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
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3, 0,
          [](const auto t) noexcept { return 3e-2 * t; }));
  problem.setPredictionPolicy(
      {.strategy =
           mfem_mgis::PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});
  // set the solver parameters
  problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                       {"AbsoluteTolerance", 1e-16},
                                       {"RelativeTolerance", 1e-16},
                                       {"MaximumNumberOfIterations", 1000}});
  problem.setSolverParameters({{"VerbosityLevel", 2},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
  //
  auto r = problem.solve(0, 1);
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
