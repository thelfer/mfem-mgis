// Benchmark MFEM Versus MFEM-MGIS (lower overhead is better)

#include <memory>
#include <cstdlib>
#include <iostream>

#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

#include <MFEMMGIS/Profiler.hxx>
#include <functional>

#include <sys/time.h>
#include <sys/resource.h>

struct TestParameters {
  const char* mesh_file = "../beam-tet.mesh";
  const char* behaviour = "Elasticity";
  const char* library = "src/libBehaviour.so";
  int order = 1;
  int refinement = 3;
  int verbosity = 0;
  int post_processing = 0;
};

void common_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 0");
  args.AddOption(&p.post_processing, "-p", "--post-processing",
                 "run post processing step");
  args.AddOption(&p.verbosity, "-v", "--verbosity", "Linear solver verbosity");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
  }
  if (mfem_mgis::getMPIrank() == 0) args.PrintOptions(std::cout);

  mfem_mgis::declareDefaultOptions(args);
}

template <typename Implementation>
void print_mesh_information(Implementation& impl) {
  using mfem_mgis::Profiler::Utils::Message;
  using mfem_mgis::Profiler::Utils::sum;
  Message("INFO: print_mesh_information");

  // getMesh
  auto mesh = impl.getFiniteElementSpace().GetMesh();

  // get the number of vertices
  int64_t numbers_of_vertices_local = mesh->GetNV();
  int64_t numbers_of_vertices = sum(numbers_of_vertices_local);

  // get the number of elements
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  // get the element size
  double h = mesh->GetElementSize(0);

  // get n dofs
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  //  int64_t unknowns = sum(unknowns_local);
  int64_t unknowns = unknowns_local;

  Message("INFO: number of vertices -> ", numbers_of_vertices);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: element size -> ", h);
  Message("INFO: Number of finite element unknowns: ", unknowns);
}

long get_memory_checkpoint() {
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  assert((test = -1) && "error: getrusage has failed");
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  return res;
};

void print_memory_footprint(std::string msg) {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6;  // conversion kb to Gb
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

int main(int argc, char* argv[]) {
  using namespace mgis::behaviour;

  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // init timers
  mfem_mgis::Profiler::timers::init_timers();

  [[maybe_unused]] double start, end;
  start = MPI_Wtime();

  // get parameters
  TestParameters p;
  double time = 0.0;
  double dt = 1.0;

  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  auto fed = mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                   {"FiniteElementFamily", "H1"},
                                   {"FiniteElementOrder", p.order},
                                   {"UnknownsSize", 3},
                                   {"NumberOfUniformRefinements", p.refinement},
                                   {"Hypothesis", "Tridimensional"},
                                   {"Parallel", true}};

  mfem_mgis::NonLinearEvolutionProblem problem(fed);
  // print_mesh_information(problem.getImplementation<true>());

  // set material properties
  problem.addBehaviourIntegrator("Mechanics", 1, p.library, p.behaviour);
  problem.addBehaviourIntegrator("Mechanics", 2, p.library, p.behaviour);

  auto& m1 = problem.getMaterial(1);
  auto& m2 = problem.getMaterial(2);

  // setting the material properties
  auto set_properties = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
    mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
  };
  set_properties(m1, 50, 50);
  set_properties(m2, 1, 1);

  // BCS

  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 1, 0));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 1, 1));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 1, 2));

  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 2, 0));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 2, 1));
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 2, 2,
          []([[maybe_unused]] const auto t) noexcept { return -1; }));

  problem.setSolverParameters({{"VerbosityLevel", 1},
                               {"RelativeTolerance", 1e-6},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 6}});

  constexpr int defaultMaxNumOfIt = 50000;  // MaximumNumberOfIterations
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(
      mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", 1e-14}});

  auto options = mfem_mgis::Parameters{{"VerbosityLevel", p.verbosity}};
  auto preconditionner =
      mfem_mgis::Parameters{{"Name", "HypreDiagScale"}, {"Options", options}};
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  // solver HyprePCG
  // problem.setLinearSolver("HypreGMRES", solverParameters);
  problem.setLinearSolver("HyprePCG", solverParameters);

  // post processing

  if (p.post_processing == 1)
    problem.addPostProcessing("ParaviewExportResults",
                              {{"OutputFileName", "Displacement"}});

  // time
  auto statistics = problem.solve(time, dt);
  if (!statistics.status) {
    mfem_mgis::Profiler::Utils::Message("INFO: FAILED");
  }
  time += dt;

  problem.update();
  if (p.post_processing == 1) problem.executePostProcessings(time, dt);

  end = MPI_Wtime();
  // printf("Duration: %1.6f s \n", end-start);fflush(stdout);
  mfem_mgis::Profiler::timers::print_timers();
  return 0;
}
