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
#if defined(MFEM_USE_MUMPS) && defined(MFEM_USE_MPI)
  auto parallel = int{1};
#else
  auto parallel = int{};
#endif
//  auto parallel = int{1};
  auto order = 1;

  // options treatment
  {
    auto args = mfem_mgis::beginParser(argc, argv);
    args->AddOption(&parallel, "-p", "--parallel",
		    "Perform parallel computations.");
    args->AddOption(&order, "-o", "--order",
		    "Finite element order (polynomial degree).");
    args->Parse();
    if (!args->Good()) {
      args->PrintUsage(std::cout);
      return EXIT_FAILURE;
    }
    args->PrintOptions(std::cout);
    mfem_mgis::endParser(args);
  }
  // loading the mesh
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
       {"Hypothesis", "Tridimensional"},
       {"Parallel", true}});

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
                                 {"RelativeTolerance", 1e-6},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10}});
    if (parallel) {
      std::cout << "MUMPS" << std::endl;
      problem.setLinearSolver("MUMPSSolver", {});
    } else {
      std::cout << "UMFPack" << std::endl;
      problem.setLinearSolver("UMFPackSolver", {});
    }
  } else {
    std::cout << "PETSc" << std::endl;
  }

  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", std::string("ssna303-displacements")}});
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}});
  
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
  return EXIT_SUCCESS;
}
