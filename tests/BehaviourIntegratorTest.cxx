/*!
 * \file   tests/BehaviourIntegratorTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/10/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem.hpp"
#include "mfem/fem/nonlinearform.hpp"
#include "MFEMMGIS/MGISIntegrator.hxx"

int main(const int argc, char** const argv) {
  const char *mesh_file = "../data/beam-tri.mesh";
  int order = 1;
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree).");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  auto mesh = std::unique_ptr<mfem::Mesh>(new mfem::Mesh(mesh_file, 1, 1));
  const auto dim = mesh->Dimension();

  if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2) {
    std::cerr << "\nInput mesh should have at least two materials and "
              << "two boundary attributes! (See schematic in ex2.cpp)\n";
    return 3;
  }
  auto fec = std::unique_ptr<mfem::H1_FECollection>(
      new mfem::H1_FECollection(order, dim));
  auto fespace =
      std::make_shared<mfem::FiniteElementSpace>(mesh.get(), fec.get(), dim);
  //
  auto i = new mfem_mgis::MGISIntegrator(fespace,
                                         mfem_mgis::Hypothesis::TRIDIMENSIONAL);
  i->addBehaviourIntegrator("SmallStrainMechanicalBehaviour", 0,
                            "src/libBehaviour.so", "Plasticity");
  //
  mfem::NonlinearForm nl(fespace.get());
  nl.AddDomainIntegrator(i);
//
  return EXIT_SUCCESS;
}