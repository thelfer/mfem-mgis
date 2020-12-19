/*!
 * \file   tests/UnilateralTensileTest.cxx
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
#include "MFEMMGIS/Material.hxx"
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
  // loading the mesh
  auto mesh = std::make_shared<mfem::Mesh>(mesh_file, 1, 1);
  if (mesh->Dimension() != dim) {
    std::cerr << "Invalid mesh dimension\n";
    return EXIT_FAILURE;
  }
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(order, dim), 3),
      mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  problem.addBehaviourIntegrator("SmallStrainMechanicalBehaviour", 1, library,
                                 behaviour);
  // materials
  auto& m1 = problem.getMaterial(1);
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
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
  //
  // Please notice there is an offset hereafter due
  // to C notation. Hence, the corresponding indices used for ess_bdr are 0,
  // 1, 4 (which corresponds to boundary attributes numbered as 1, 2, 5). This
  // means that we want to fix dirichlet boundary conditions on "xz0", "xy0"
  // and "yz0".
  auto fixed_dirichlet_dofs = mfem::Array<mfem_mgis::size_type>{};
  auto append_dof = [&](const auto boundary, const auto normal) {
    auto tmp = mfem::Array<mfem_mgis::size_type>{};
    auto boundaries_markers =
        mfem::Array<mfem_mgis::size_type>(mesh->bdr_attributes.Max());
    boundaries_markers = 0;
    boundaries_markers[boundary] = 1;
    problem.getFiniteElementSpace().GetEssentialTrueDofs(boundaries_markers,
                                                         tmp, normal);
    fixed_dirichlet_dofs.Append(tmp);
    return tmp;
  };
  append_dof(0, 1);                     // xz0
  append_dof(1, 2);                     // xy0
  append_dof(4, 0);                     // yz0
  auto yz1_ux_dofs = append_dof(2, 0);  // yz1
  problem.SetEssentialTrueDofs(fixed_dirichlet_dofs);
  // solving the problem
  mfem::UMFPackSolver lsolver;
  auto& solver = problem.getSolver();
  solver.iterative_mode = true;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(0);
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(0);
  solver.SetMaxIter(10);
  // vtk export
  mfem::ParaViewDataCollection paraview_dc("UnilateralTensileTestOutput",
                                           mesh.get());
  const auto vo =
      mgis::behaviour::getVariableOffset(m1.b.isvs, isv_name, m1.b.hypothesis);
  auto exx = std::vector<mfem_mgis::real>{};
  auto eyy = std::vector<mfem_mgis::real>{};
  auto sxx = std::vector<mfem_mgis::real>{};
  auto v = std::vector<mfem_mgis::real>{};
  // loop over time step
  exx.push_back(m1.s0.gradients[0]);
  eyy.push_back(m1.s0.gradients[1]);
  sxx.push_back(m1.s0.thermodynamic_forces[0]);
  v.push_back(m1.s0.internal_state_variables[vo]);
  const auto nsteps = mfem_mgis::size_type{100};
  const auto dt = mfem_mgis::real{1} / nsteps;
  const auto u_t = [](const auto t) {
    if (t < 0.3) {
      return 3e-2 * t;
    } else if (t < 0.6) {
      return 0.009 - 0.1 * (t - 0.3);
    }
    return -0.021 + 0.1 * (t - 0.6);
  };
  auto t = mfem_mgis::real{0};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    // setting the boundary values
    auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
    u1.Print();
    const auto u = u_t(t + dt);
    for (mfem_mgis::size_type idx = 0; idx != yz1_ux_dofs.Size(); ++idx) {
      u1[yz1_ux_dofs[idx]] = u;
    }
    // resolution
    problem.solve(dt);
    problem.update();
    t += dt;
    //
    exx.push_back(m1.s1.gradients[0]);
    eyy.push_back(m1.s1.gradients[1]);
    sxx.push_back(m1.s1.thermodynamic_forces[0]);
    v.push_back(m1.s1.internal_state_variables[vo]);
    // recover the solution as a grid function
    mfem::GridFunction x(&problem.getFiniteElementSpace());
    x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
    x.SetFromTrueVector();
    paraview_dc.RegisterField("u", &x);
    paraview_dc.SetCycle(i);
    paraview_dc.SetTime(t);
    paraview_dc.Save();
  }
  // save the traction curve
  std::ofstream out("UnilateralTensileTest-" + std::string(behaviour) + ".txt");
  out.precision(14);
  for (std::vector<mfem_mgis::real>::size_type i = 0; i != exx.size(); ++i) {
    out << exx[i] << " " << eyy[i] << " " << sxx[i] << " " << v[i] << '\n';
  }
  // comparison to reference results
  constexpr const auto eps = mfem_mgis::real(1.e-10);
  constexpr const auto E = mfem_mgis::real(70.e9);
  auto success = true;
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
  std::ifstream in(reference_file);
  for (std::vector<mfem_mgis::real>::size_type i = 0; i != exx.size(); ++i) {
    auto exx_ref = mfem_mgis::real{};
    auto eyy_ref = mfem_mgis::real{};
    auto sxx_ref = mfem_mgis::real{};
    auto v_ref = mfem_mgis::real{};
    in >> exx_ref >> eyy_ref >> sxx_ref >> v_ref;
    check(sxx[i], sxx_ref, E * eps, "invalid stress value");
    check(eyy[i], eyy_ref, eps, "invalid transverse strain");
    check(v[i], v_ref, eps, "invalid internal state variable");
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
