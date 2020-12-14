/*!
 * \file   tests/PeriodicTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/10/2020
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
  auto tcase = 1;
  auto static_cond = false;
  auto visualization = 1;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library", "Material library.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&tcase, "-t", "--tcase",
                 "identifier of the case : Exx->1, Eyy->2, Ezz->3, Exy->4, "
                 "Eyz->5, Exz->6");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  args.PrintOptions(std::cout);
  // parsing the mesh
  auto mesh = std::make_shared<mfem::Mesh>(mesh_file, 1, 1);
  const auto dim = mesh->Dimension();
  std::cout << "dim: " << dim << std::endl;
  // creating the finite element space
  auto const fec = std::make_shared<mfem::H1_FECollection>(order, dim);
  auto fespace =
      std::make_shared<mfem::FiniteElementSpace>(mesh.get(), fec.get(), dim);
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem p(
      fespace, mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  p.addBehaviourIntegrator("SmallStrainMechanicalBehaviour", 1, library,
                           "Elasticity");
  p.addBehaviourIntegrator("SmallStrainMechanicalBehaviour", 2, library,
                           "Elasticity");
  // materials
  auto& m1 = p.getMaterial(1);
  auto& m2 = p.getMaterial(2);
  // setting the material properties
  auto set_properties = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
  };
  set_properties(m1, 100, 75);
  set_properties(m2, 200, 150);
  //
  auto set_temperature = [](auto& m) {
    mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
  };
  set_temperature(m1);
  set_temperature(m2);
  // macroscopic strain
  std::vector<double> e = {1, 0, 0, 0, 0, 0};
  m1.setMacroscopicGradients(e);
  m2.setMacroscopicGradients(e);
  // Impose no displacement on the first node
  // which needs to be on x=xmin or x=xmax axis.
  // ux=0, uy=0, uz=0 on this point.
  int ndof = fespace->GetTrueVSize() / dim;
  mfem::Array<int> ess_tdof_list;
  ess_tdof_list.SetSize(dim);
  for (int k = 0; k < dim; k++) {
    int tgdof = 0 + k * ndof;
    ess_tdof_list[k] = tgdof;
  }
  p.SetEssentialTrueDofs(ess_tdof_list);
  // solving the problem
  mfem::UMFPackSolver lsolver;
  auto& solver = p.getSolver();
  solver.iterative_mode = false;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(1);  // print Newton iterations
  solver.SetRelTol(1e-4);
  solver.SetAbsTol(1e-10);
  solver.SetMaxIter(10);
  p.solve(1);
  //
  auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
  mfem::GridFunction u1_gf(fespace.get());
  u1_gf.MakeTRef(fespace.get(), u1, 0);
  u1_gf.SetFromTrueVector();
  mfem::ParaViewDataCollection paraview_dc("per", mesh.get());
  paraview_dc.RegisterField("u", &u1_gf);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);
  paraview_dc.Save();
  return EXIT_SUCCESS;
}
