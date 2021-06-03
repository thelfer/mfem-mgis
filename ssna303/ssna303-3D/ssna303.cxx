/*!
 * \file   ssna303.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/linalg/hypre.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"

#define PRINT_DEBUG (std::cout <<  __FILE__ << ":" <<  __LINE__ << std::endl)

int main(int argc, char** argv) {
	
  mfem_mgis::initialize(argc, argv);
  bool parallel = true;
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  auto solver = "HyprePCG";
  auto preconditioner = "";
  auto ref_para = 0;
  auto ref_seq = 0;
  auto order = 1;
  // initialization for the MPI 
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD,&rank);


  //file creation 
  std::string const myFile("/home/hc265945/spack_codes/mfem-mgis/ssna303/ssna303-3D/Test_Ssna303/Grandes_Perturbations/Without_Preconditioner/HyprePCG/HyprePCG.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&solver,"-s", "--solver",
                 "Solver of the Problem.");
  args.AddOption(&preconditioner,"-p", "--preconditioner",
                 "Preconditioner uses for the Problem.");
  args.AddOption(&ref_para,"-rp", "--refinement_parallel",
                 "Number of Refinement for parallel call.");
  args.AddOption(&ref_seq,"-rs", "--refinement_sequential",
                 "Number of Refinement for sequential call.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);

  // loading the mesh
  {
  const auto main_timer = mfem_mgis::getTimer("main_timer");
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
       {"NumberOfUniformRefinements", parallel ? ref_para : ref_seq},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", true}});

  auto mesh = problem.getImplementation<true>().getFiniteElementSpace().GetMesh();
  //get the number of vertices
  int numbers_of_vertices = mesh->GetNV();
  //get the number of elements
  int numbers_of_elements = mesh->GetNE();
  //get the element size
  double h = mesh->GetElementSize(0);

  // 2 1 "Volume"
  problem.addBehaviourIntegrator("Mechanics", 1, library, behaviour);
  // materials
  auto& m1 = problem.getMaterial(1);
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
  // boundary conditions

  // 3 LowerBoundary
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3, 1));
  // 4 SymmetryPlane1
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 4, 0));
  // 5 SymmetryPlane2
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 5, 2));
  // 2 UpperBoundary 
  problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 2, 1,
          [](const auto t) {
            const auto u = 6e-3 * t;
            return u;
          }));

  // solving the problem
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-6},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});

  // selection of the linear solver without preconditioner
  if (solver == ""){
    return EXIT_FAILURE;
  }
  if ( preconditioner == ""){
  problem.setLinearSolver(solver,  {{"VerbosityLevel", 0},
                                   {"RelativeTolerance", 1e-12},
                                   {"MaximumNumberOfIterations",300}});
  }
  else{
  // with the HypreBoomerAMG preconditioner
//  auto prec_none = mfem_mgis::Parameters{{"Name", "None"}};
  auto prec_boomer =
      mfem_mgis::Parameters{{"Name", preconditioner},  //
                           {"Options", mfem_mgis::Parameters{
//                           {"Strategy", "Elasticity"},
                           {"VerbosityLevel", 0}}}};

   problem.setLinearSolver(solver, {{"VerbosityLevel", 0},
                                          {"AbsoluteTolerance", 1e-12},
                                          {"RelativeTolerance", 1e-12},
                                          {"MaximumNumberOfIterations", 300},
                                          {"Preconditioner", prec_boomer}});
  }
	// print on file
  out << "SetLinearSolver" << std::endl;
  out << "VerbosityLevel = " << 0 << std::endl;
  out << "RelativeTolerance = " << 1e-12 << std::endl;
  out << "MaximumNumberOfIterations = " << 300 << std::endl;
  out << "Preconditioner = " << preconditioner << std::endl;
  out << "taille_maille = " << h << std::endl;
  out << "1/h = " << 1/h << std::endl;
  out << " nbr_ref_parallel = " << ref_para << std::endl;
  out << " nbr_ref_sequential = " << ref_seq << std::endl;
  out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
  out << " numbers_of_elements = " << numbers_of_elements << std::endl;

  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", std::string("ssna303-displacements-HPCG_WP_0")}});
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force_HPCG_WP_0.txt"}});
  
  // loop over time step
  const auto nsteps = mfem_mgis::size_type{50};
  const auto dt = mfem_mgis::real{1} / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
    // resolution
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto niter  = mfem_mgis::size_type{0};
    while (nsteps != 0) {
      bool converged = true;
      try {
        problem.solve(ct, dt2);
      } catch (std::runtime_error&) {
        converged = false;
      }
      if (converged) {
        --nsteps;
        ct += dt2;
      } else {
        std::cout << "\nsubstep: " << niter << '\n';
        nsteps *= 2;
        dt2 /= 2;
        ++niter;
        problem.revert();
        if (niter == 10) {
          mgis::raise("maximum number of substeps");
        }
      }
    }
    problem.executePostProcessings(t, dt);
    problem.update();
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  }
  mfem_mgis::Profiler::getProfiler().print(out);
  return EXIT_SUCCESS;
  MPI_Finalize();
}
