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
  // setting the material properties
  auto set_properties = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
  };
  set_properties(p.getMaterial(1), 100, 75);
  set_properties(p.getMaterial(2), 100, 75);
  //
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
  p.SetEssentialBC(ess_tdof_list);
  // solving the problem
  p.solve(1);
  return EXIT_SUCCESS;
}
