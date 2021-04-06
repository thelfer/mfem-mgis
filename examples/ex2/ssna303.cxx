/*!
 * \file   ssna3030.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

void exit_on_failure() {
  MPI_Finalize();
  std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
   // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  constexpr const auto dim = mfem_mgis::size_type{2};
  const char* mesh_file = "ssna303.msh";
  const char* behaviour = "Plasticity";
  const char* library = "src/libBehaviour.so";
  auto order = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    exit_on_failure();
  }
  args.PrintOptions(std::cout);
  // loading the mesh
  auto smesh = std::make_shared<mfem::Mesh>(mesh_file, 1, 1);
  auto mesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *smesh); //
  //  auto mesh = smesh; // Uncomment this line for serial version
  if (mesh->Dimension() != dim) {
    std::cerr << "Invalid mesh dimension\n";
    exit_on_failure();
  }
  //
  const auto& bdr_attributes = mesh->bdr_attributes;
  std::cout << "boundaries: ";
  for (mfem_mgis::size_type i = 0; i != bdr_attributes.Size(); ++i) {
    std::cout << " " << bdr_attributes[i];
  }
  std::cout << '\n';
  //
  const auto& material_attributes = mesh->attributes;
  std::cout << "materials: ";
  for (mfem_mgis::size_type i = 0; i != material_attributes.Size(); ++i) {
    std::cout << " " << material_attributes[i];
  }
  std::cout << '\n';
  // building the non linear problem
  {
    mfem_mgis::NonLinearEvolutionProblem problem(
        std::make_shared<mfem_mgis::FiniteElementDiscretization>(
            mesh, std::make_shared<mfem::H1_FECollection>(order, dim), dim),
        mgis::behaviour::Hypothesis::PLANESTRAIN);
    // 2 1 "Volume"
    problem.addBehaviourIntegrator("Mechanics", 1, library, behaviour);
    // materials
    auto& m1 = problem.getMaterial(1);
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
    // boundary conditions
    // 1 3 "LowerBoundary"
    problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(), 3, 1));
    // 1 4 "SymmetryAxis "
    problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(), 4, 0));
    // 1 2 "UpperBoundary"
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
  // selection of the linear solver
#if defined (MFEM_USE_MUMPS)
    problem.setLinearSolver("MUMPSSolver", {});
#else
    problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                         {"RelativeTolerance", 1e-12},
                                         {"MaximumNumberOfIterations", 300}});
#endif
    // vtk export
    problem.addPostProcessing("ParaviewExportResults",
                              {{"OutputFileName", std::string("ssna303-displacements")}});
//    problem.addPostProcessing("ComputeResultantForceOnBoundary",
//                              {{"Boundary", 2}, {"OutputFileName", "force.txt"}});
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
        auto converged = true;
        try {
          problem.solve(ct, dt2);
        } catch (std::runtime_error&) {
          converged = false;
        }
        //      const auto converged = problem.solve(ct, dt2);
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
  MPI_Finalize();
  return EXIT_SUCCESS;
}
