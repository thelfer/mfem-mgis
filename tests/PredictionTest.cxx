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

template <bool parallel>
void export_prediction(mfem_mgis::NonLinearEvolutionProblem &p,
                       mfem_mgis::GridFunction<parallel> &du) {
  auto &fed = p.getFiniteElementDiscretization();
  auto exporter = mfem::ParaViewDataCollection{
      parallel ? "result-prediction-test-parallel" : "result-prediction-test"};
  exporter.SetMesh(&(fed.getMesh<parallel>()));
  exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  exporter.RegisterField("DisplacementPrediction", &du);
  exporter.SetCycle(1);
  exporter.SetTime(1);
  exporter.Save();
}

template <bool parallel>
[[nodiscard]] mfem_mgis::GridFunction<parallel> computePrediction(
    mfem_mgis::Context &ctx, mfem_mgis::NonLinearEvolutionProblem &p) {
  auto &fed = p.getFiniteElementDiscretization();
  auto &fespace = fed.getFiniteElementSpace<parallel>();
  const auto &u0 = p.getUnknowns(mfem_mgis::bts);
  p.setup(0, 1);
  auto success =
      p.integrate(u0, mfem_mgis::IntegrationType::PREDICTION_ELASTIC_OPERATOR);
  if constexpr (parallel) {
    MPI_Allreduce(MPI_IN_PLACE, &success, 1, MPI_C_BOOL, MPI_LAND,
                  fespace.GetComm());
  }
  if (!success) {
    mfem_mgis::getOutputStream() << "integration failure\n";
    std::exit(EXIT_FAILURE);
  }
  auto operators = p.getLinearizedOperators(ctx, u0);
  if (isInvalid(operators)) {
    mfem_mgis::getOutputStream() << ctx.getErrorMessage() << '\n';
    std::exit(EXIT_FAILURE);
  }
  //
  auto essential_dofs = p.getEssentialDegreesOfFreedom();
  auto edofs_list = mfem::Array<mfem_mgis::size_type>(essential_dofs.data(),
                                                      essential_dofs.size());
  mfem_mgis::GridFunction<parallel> du(&fespace);
  if constexpr (parallel) {
    auto du_values = mfem::Vector(fespace.GetTrueVSize());
    du_values = 0.0;
    //
    for (const auto &bc : p.getDirichletBoundaryConditions()) {
      bc->setImposedValuesIncrements(du_values, 0, 1);
    }
    du.Distribute(du_values);
  } else {
    for (const auto &bc : p.getDirichletBoundaryConditions()) {
      bc->setImposedValuesIncrements(du, 0, 1);
    }
  }
  //
  auto a = mfem_mgis::BilinearForm<parallel>(&fespace);
  a.AddDomainIntegrator(operators->K.release());
  a.Assemble();
  auto b = mfem_mgis::LinearForm<parallel>(&fespace);
  b.AddDomainIntegrator(operators->Fi.release());
  b.Assemble();
  //
  auto A = [] {
    if constexpr (parallel) {
      return mfem::HypreParMatrix{};
    } else {
      return mfem::SparseMatrix{};
    }
  }();
  mfem::Vector B, X;
  //
  a.FormLinearSystem(edofs_list, du, b, A, X, B);
  //
  auto ls = [&fespace, &ctx] {
    auto &f = mfem_mgis::LinearSolverFactory<parallel>::getFactory();
    if constexpr (parallel) {
      return f.generate(ctx, "CGSolver", fespace,
                        {{"VerbosityLevel", 1},
                         {"AbsoluteTolerance", 1e-12},
                         {"RelativeTolerance", 1e-12},
                         {"MaximumNumberOfIterations", 5000}});
    } else {
      return f.generate(ctx, "UMFPackSolver", fespace, {});
    }
  }();
  if (isInvalid(ls)) {
    mfem_mgis::getOutputStream() << ctx.getErrorMessage() << '\n';
    std::exit(EXIT_FAILURE);
  }
  ls.linear_solver->SetOperator(A);
  ls.linear_solver->Mult(B, X);
  a.RecoverFEMSolution(X, b, du);
  return du;
}  // end of computePrediction

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
         {"NumberOfUniformRefinements", 2},  // faster for testing
         //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
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
    //
    if (parallel) {
#ifdef MFEM_USE_MPI
      auto du = computePrediction<true>(ctx, problem);
      export_prediction<true>(problem, du);
#else  /* MFEM_USE_MPI */
      return EXIT_FAILURE;
#endif /* MFEM_USE_MPI */
  } else {
    auto du = computePrediction<false>(ctx, problem);
    export_prediction<false>(problem, du);
  }
  return EXIT_SUCCESS;
}
