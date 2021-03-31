/*!
 * \file   tests/UniaxialTensileTest.cxx
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
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

int main(const int argc, char** const argv) {
  constexpr const auto dim = mfem_mgis::size_type{3};
  const char* mesh_file = nullptr;
  const char* behaviour = nullptr;
  const char* library = nullptr;
  const char* reference_file = nullptr;
  const char* isv_name = nullptr;
  auto order = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&reference_file, "-r", "--reference-file", "Reference file.");
  args.AddOption(&behaviour, "-b", "--behaviour", "Name of the behaviour.");
  args.AddOption(&isv_name, "-v", "--internal-state-variable",
                 "Internal variable name to be post-processed.");
  args.AddOption(&library, "-l", "--library", "Material library.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.Parse();
  if ((!args.Good()) || (mesh_file == nullptr) || (reference_file == nullptr) ||
      (isv_name == nullptr) || (behaviour == nullptr)) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"Hypothesis", "Tridimensional"},
       {"UnknownsSize", dim},
       {"Parallel", false}});
  problem.addBehaviourIntegrator("Mechanics", 1, library, behaviour);
  // materials
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
  // solving the problem
  problem.setSolverParameters({{"VerbosityLevel", 0},
                               {"RelativeTolerance", 1e-12},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 10}});
#ifdef MFEM_USE_SUITESPARSE
  problem.setLinearSolver("UMFPackSolver", {});
#else
  problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
                                       {"RelativeTolerance", 1e-12},
                                       {"MaximumNumberOfIterations", 300}});
#endif
  // vtk export
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "UniaxialTensileTestOutput-" +
                                                    std::string(behaviour)}});
  //
  const auto vo =
      mgis::behaviour::getVariableOffset(m1.b.isvs, isv_name, m1.b.hypothesis);
  auto g0 = std::vector<mfem_mgis::real>{};
  auto g1 = std::vector<mfem_mgis::real>{};
  auto tf0 = std::vector<mfem_mgis::real>{};
  auto v = std::vector<mfem_mgis::real>{};
  // loop over time step
  g0.push_back(m1.s0.gradients[0]);
  g1.push_back(m1.s0.gradients[1]);
  tf0.push_back(m1.s0.thermodynamic_forces[0]);
  v.push_back(m1.s0.internal_state_variables[vo]);
  const auto nsteps = mfem_mgis::size_type{100};
  const auto dt = mfem_mgis::real{1} / nsteps;
  auto t = mfem_mgis::real{0};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    // resolution
    if (!problem.solve(t, dt)) {
      return EXIT_FAILURE;
    }
    problem.executePostProcessings(t, dt);
    problem.update();
    t += dt;
    //
    g0.push_back(m1.s1.gradients[0]);
    g1.push_back(m1.s1.gradients[1]);
    tf0.push_back(m1.s1.thermodynamic_forces[0]);
    v.push_back(m1.s1.internal_state_variables[vo]);
  }
  // save the traction curve
  std::ofstream out("UniaxialTensileTest-" + std::string(behaviour) + ".txt");
  out.precision(14);
  for (std::vector<mfem_mgis::real>::size_type i = 0; i != g0.size(); ++i) {
    out << g0[i] << " " << g1[i] << " " << tf0[i] << " " << v[i] << '\n';
  }
  // comparison to reference results
  auto success = true;
  std::ifstream in(reference_file);
  if (in) {
    constexpr const auto eps = mfem_mgis::real(1.e-10);
    constexpr const auto E = mfem_mgis::real(70.e9);
    auto check = [&success](const auto cv,  // computed value
                            const auto rv,  // reference value,
                            const auto ev, const auto msg) {
      const auto e = std::abs(cv - rv);
      if (e > ev) {
        std::cerr << "test failed (" << msg << ", " << cv << " vs " << rv
                  << ", error " << e << ")\n";
        success = false;
      }
    };
    for (std::vector<mfem_mgis::real>::size_type i = 0; i != g0.size(); ++i) {
      auto g0_ref = mfem_mgis::real{};
      auto g1_ref = mfem_mgis::real{};
      auto tf0_ref = mfem_mgis::real{};
      auto v_ref = mfem_mgis::real{};
      in >> g0_ref >> g1_ref >> tf0_ref >> v_ref;
      check(tf0[i], tf0_ref, E * eps, "invalid stress value");
      check(g1[i], g1_ref, eps, "invalid transverse strain");
      check(v[i], v_ref, eps, "invalid internal state variable");
    }
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
