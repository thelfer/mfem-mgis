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
  const char* behaviour = "Elasticity";
  const char* library = "src/libBehaviour.so";
  const char* petscrc_file = "petsc_file";
  auto parallel = int{1};
  auto order = 1;
  auto refinement =0;
  
  //file creation 
  std::string const myFile("./data.txt");
  std::ofstream out(myFile.c_str());

  // options treatment
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);// PETSc Initialize 
  if (!mfem_mgis::usePETSc())  
    mfem_mgis::setPETSc(petscrc_file);
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

  {
  // loading the mesh and timer
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
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

  // print on file
  out << " USE_PETSc = " << mfem_mgis::usePETSc() << std::endl;
  out << " taille_maille = " << h << std::endl;
  out << " 1/h = " << 1/h << std::endl;
  out << " nbr_refinement = " << refinement << std::endl;
  out << " numbers_of_vertices = " << numbers_of_vertices << std::endl;
  out << " numbers_of_elements = " << numbers_of_elements << std::endl;
  
  // loop over time step
//  const auto nsteps = mfem_mgis::size_type{100};
//  const auto dt = mfem_mgis::real{1} / nsteps;
  const auto nsteps = mfem_mgis::size_type{2};
  const auto dt = mfem_mgis::real{0.001};
  auto t = mfem_mgis::real{0};
  bool converged = problem.solve(t, dt);
  if (!converged) { mgis::raise("error, this problem should converge"); }
  problem.update();
  problem.executePostProcessings(t, dt);
  mfem_mgis::finalize();
  return EXIT_SUCCESS;
}
