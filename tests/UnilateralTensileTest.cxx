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
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  auto order = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library", "Material library.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);
  // parsing the mesh
  auto mesh = std::make_shared<mfem::Mesh>(mesh_file, 1, 1);
  const auto dim = mesh->Dimension();
  // creating the finite element space
  auto const fec = std::make_shared<mfem::H1_FECollection>(order, dim);
  auto fespace =
      std::make_shared<mfem::FiniteElementSpace>(mesh.get(), fec.get(), dim);
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem p(
      fespace, mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  p.addBehaviourIntegrator("SmallStrainMechanicalBehaviour", 1, library,
                           "Plasticity");
  // materials
  auto& m1 = p.getMaterial(1);
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
    auto tmp = mfem::Array<int>{};
    auto boundaries_markers = mfem::Array<int>(mesh->bdr_attributes.Max());
    boundaries_markers = 0;
    boundaries_markers[boundary] = 1;
    fespace->GetEssentialTrueDofs(boundaries_markers, tmp, normal);
    fixed_dirichlet_dofs.Append(tmp);
    return tmp;
  };
  append_dof(0, 1);  // xz0
  append_dof(1, 2);  // xy0
  append_dof(4, 0);  // yz0
  auto yz1_ux_dofs = append_dof(2, 0);  // yz1
  p.SetEssentialTrueDofs(fixed_dirichlet_dofs);
  // solving the problem
  mfem::UMFPackSolver lsolver;
  auto& solver = p.getSolver();
  solver.iterative_mode = false;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(1);  // print Newton iterations
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(0);
  solver.SetMaxIter(10);
  // vtk export
  mfem::ParaViewDataCollection paraview_dc("UnilateralTensileTestOutput",
                                           mesh.get());
  // loop over time step
  auto t = mfem_mgis::real{0};
  auto dt = mfem_mgis::real{1};
  for (mfem_mgis::size_type i = 0; i != 10; ++i) {
    // setting the boundary values
    auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    const auto u = 1e-2 * (t + dt) / 10;
    for (mfem_mgis::size_type idx = 0; idx != yz1_ux_dofs.Size(); ++idx) {
      u1[yz1_ux_dofs[idx]] = u;
    }
    // resolution
    p.solve(dt);
    p.update();
    t += dt;
    // recover the solution as a grid function
    mfem::GridFunction x(fespace.get());
    x.MakeTRef(fespace.get(), u1, 0);
    x.SetFromTrueVector();
    paraview_dc.RegisterField("u", &x);
    paraview_dc.SetCycle(i);
    paraview_dc.SetTime(t);
    paraview_dc.Save();
  }
  return EXIT_SUCCESS;
}
