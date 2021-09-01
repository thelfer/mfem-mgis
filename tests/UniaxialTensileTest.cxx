/*!
 * \file   tests/UniaxialTensileTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
#ifdef DO_USE_MPI
  static constexpr const auto parallel = true;
#else
  static constexpr const auto parallel = false;
#endif
  constexpr const auto dim = mfem_mgis::size_type{3};
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  auto success = true;
  {
    const auto main_timer = mfem_mgis::getTimer("main");
    // building the non linear problem
    mfem_mgis::NonLinearEvolutionProblem problem(
        {{"MeshFileName", parameters.mesh_file},
         {"FiniteElementFamily", "H1"},
         {"FiniteElementOrder", parameters.order},
         {"UnknownsSize", dim},
         {"NumberOfUniformRefinements", parallel ? 2 : 0},
         {"Hypothesis", "Tridimensional"},
         {"Parallel", parallel}});
    // materials
    problem.addBehaviourIntegrator("Mechanics", 1, parameters.library,
                                   parameters.behaviour);
    auto& m1 = problem.getMaterial(1);
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
    if (m1.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
      std::array<mfem_mgis::real, 9u> r = {0, 1, 0,  //
                                           1, 0, 0,  //
                                           0, 0, 1};
      m1.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    }
    // boundary conditions
    // Determine the list of true (i.e. parallel conforming) essential
    // boundary dofs. In this example, the boundary conditions are defined by
    // marking only boundary attribute 1 (xz0) 2 (xy0) and 5 (yz0) from the
    // mesh as essential and converting it to a list of true dofs.
    // In the mesh file, we have both an index and a name for physical entities
    // :
    //    $PhysicalNames
    //    7
    //    2 1 "xz0"
    //    2 2 "xy0"
    //    2 3 "yz1"
    //    2 4 "xz1"
    //    2 5 "yz0"
    //    2 6 "xy1"
    //    3 7 "mat"
    //    $EndPhysicalNames
    // Only the index is used in this C++ code for manipulating related dof.
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
            [](const auto t) {
              if (t < 0.3) {
                return 3e-2 * t;
              } else if (t < 0.6) {
                return 0.009 - 0.1 * (t - 0.3);
              }
              return -0.021 + 0.1 * (t - 0.6);
            }));
    // set the solver parameters
    mfem_mgis::unit_tests::setLinearSolver(problem, parameters);
    problem.setSolverParameters({{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-12},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10}});
    // vtk export
    problem.addPostProcessing(
        "ParaviewExportResults",
        {{"OutputFileName",
          "UniaxialTensileTestOutput-" + std::string(parameters.behaviour)}});
    const auto& b = problem.getMaterial(1).b;
    if((b.btype==mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR)&&
       (b.kinematic==mgis::behaviour::Behaviour::SMALLSTRAINKINEMATIC)){
      problem.addPostProcessing(
				"ParaviewExportIntegrationPointResultsAtNodes",
				{{"OutputFileName", "UniaxialTensileTestIntegrationPointOutput-" +
				  std::string(parameters.behaviour)},
				 {"Materials", {1}},
				 {"Results", {"Strain"}}});
    }
    // solving the problem in 100 time steps
    auto r = mfem_mgis::unit_tests::solve(problem, parameters, 0, 1, 100);
    // save the results curve
    mfem_mgis::unit_tests::saveResults(
        "UniaxialTensileTest-" + std::string(parameters.behaviour) + ".txt", r);
    // compare to reference files
    constexpr const auto eps = mfem_mgis::real(1.e-10);
    constexpr const auto E = mfem_mgis::real(70.e9);
    success =
        mfem_mgis::unit_tests::checkResults(r, m1, parameters, eps, E * eps);
  }
  //  mfem_mgis::Profiler::getProfiler().print(mfem_mgis::getOutputStream());
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
