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
  constexpr const auto eps = mfem_mgis::real{1e-10};
  constexpr const auto xthr = mfem_mgis::real(1) / 2;
  std::array<void (*)(const mfem::Vector&, mfem::Vector&), 6u> solutions = {
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector&, mfem::Vector& u) { u = mfem_mgis::real{}; }};
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  auto order = 1;
  auto tcase = 0;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&library, "-l", "--library", "Material library.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
  args.PrintOptions(std::cout);
  if ((tcase < 0) || (tcase > 5)) {
    std::cerr << "Invalid test case\n";
    return EXIT_FAILURE;
  }
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
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
  if (tcase < 3) {
    e[tcase] = 1;
  } else {
    e[tcase] = 1.41421356237309504880 / 2;
  }
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
  solver.SetRelTol(1e-12);
  solver.SetAbsTol(1e-12);
  solver.SetMaxIter(10);
  p.solve(1);
  // recover the solution as a grid function
  auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
  mfem::GridFunction x(fespace.get());
  x.MakeTRef(fespace.get(), u1, 0);
  x.SetFromTrueVector();
  // comparison to analytical solution
  mfem::VectorFunctionCoefficient sol_coef(dim, solutions[tcase]);
  const auto error = x.ComputeL2Error(sol_coef);
  if (error > eps) {
    std::cerr << "Error is greater than threshold (" << error << " > " << eps << ")\n";
    return EXIT_FAILURE;
  } else {
    std::cerr << "Error is lower than threshold (" << error << " < " << eps << ")\n";
  }
  // exporting the results
  mfem::ParaViewDataCollection paraview_dc(
      "PeriodicTest-" + std::to_string(tcase), mesh.get());
  paraview_dc.RegisterField("u", &x);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);
  paraview_dc.Save();
  return EXIT_SUCCESS;
}
