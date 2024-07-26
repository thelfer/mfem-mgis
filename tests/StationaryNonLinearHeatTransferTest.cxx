/*!
 * \file   tests/StationaryNonLinearHeatTransferTest.cxx
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
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "UnitTestingUtilities.hxx"

int main(int argc, char** argv) {
  auto parameters = mfem_mgis::unit_tests::TestParameters{};
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::unit_tests::parseCommandLineOptions(parameters, argc, argv);
  if (parameters.isv_name != nullptr) {
    mfem_mgis::abort("no internal state variable expected");
  }
  auto success = true;
  {
    // building the non linear problem
    mfem_mgis::NonLinearEvolutionProblem problem(
        {{"MeshFileName", parameters.mesh_file},
         {"FiniteElementFamily", "H1"},
         {"FiniteElementOrder", parameters.order},
         {"UnknownsSize", 1},
         {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
         {"Hypothesis", "Tridimensional"},
         {"Parallel", bool(parameters.parallel)}});
    //
    problem.getUnknownsAtBeginningOfTheTimeStep() = 293.15;
    problem.getUnknownsAtEndOfTheTimeStep() = 293.15;
    // materials
    problem.addBehaviourIntegrator("StationaryNonLinearHeatTransfer", 1,
                                   parameters.library, parameters.behaviour);
    auto& m1 = problem.getMaterial(1);
    auto T = std::vector<mfem_mgis::real>(m1.n, 293.15);
    mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", T);
    mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", T);
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
            problem.getFiniteElementDiscretizationPointer(), 5, 0,
            [](const auto) { return 293.15; }));
    problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(), 3, 0,
            [](const auto t) { return 293.15 + (893.15 - 293.15) * t; }));
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
    // solving the problem in 100 time steps
    auto r = mfem_mgis::unit_tests::solve(problem, parameters, 0, 1, 100);
    // save the results curve
    mfem_mgis::unit_tests::saveResults(
        "UniaxialTensileTest-" + std::string(parameters.behaviour) + ".txt", r);
    //     // compare to reference files
    //     constexpr const auto eps = mfem_mgis::real(1.e-10);
    //     constexpr const auto E = mfem_mgis::real(70.e9);
    //     success =
    //         mfem_mgis::unit_tests::checkResults(r, m1, parameters, eps, E *
    //         eps);
  }
  mfem_mgis::Profiler::timers::print_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
