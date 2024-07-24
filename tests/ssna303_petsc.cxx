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

#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#ifdef MFEM_USE_MUMPS
#include "mfem/linalg/mumps.hpp"
#endif /* MFEM_USE_MUMPS */

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
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = "ssna303_3d.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  const char* petscrc_file = "rc_ex10p";
  int parallel = int{1};
  auto order = 1;
  auto refinement =0;
  
  //file creation 
  std::string const myFile("./data.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);// PETSc Initialize 
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library", "behaviour law lib.");
  args.AddOption(&petscrc_file, "-pf", "--petsc-file",
                 "path to the petscfile.");
  args.AddOption(&parallel, "-p", "--parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&refinement,"-r", "--refinement",
                 "Number of Refinement.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);

  if (!mfem_mgis::usePETSc())  
    mfem_mgis::setPETSc(petscrc_file);

  {
  // loading the mesh and timer
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", bool(parallel)}});

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

  // solving the problem without petsc
  if (!mfem_mgis::usePETSc()) {
    problem.setSolverParameters({{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-8},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 20}});
    if (parallel) {
      std::cout << "MUMPS" << std::endl;
      problem.setLinearSolver("MUMPSSolver", {});
    } else {
        std::cout << "UMFSolver" << std::endl;
        problem.setLinearSolver("UMFPackSolver", {});
    }
  }

  // print on file
  out << " USE_PETSc = " << mfem_mgis::usePETSc() << std::endl;

  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", std::string("ssna303-displacements")}});
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}});
  
  // loop over time step
//  const auto nsteps = mfem_mgis::size_type{100};
//  const auto dt = mfem_mgis::real{1} / nsteps;
  const auto nsteps = mfem_mgis::size_type{2};
  const auto dt = mfem_mgis::real{0.001};
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
      bool converged = problem.solve(ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        problem.update();
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
    t += dt;
    ++iteration;
    std::cout << '\n';
  }
  }
  mfem_mgis::finalize();
  return EXIT_SUCCESS;
}
