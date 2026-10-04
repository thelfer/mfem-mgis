/*!
 * \file   HypreBoomerAMGTest.cxx
 * \brief  Test of the system BoomerAMG preconditioner built by the linear
 * solver factory
 * \date   02/10/2026
 */

#include <cstdlib>
#include <iostream>
#include "mfem/mesh/pmesh.hpp"
#include "mfem/fem/pbilinearform.hpp"
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pgridfunc.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

/*!
 * \return the number of iterations of a preconditioned conjugate gradient
 * solving twice the given system, as two Newton iterations would do
 * \param[in] fespace: finite element space
 * \param[in] A: matrix
 * \param[in] B: right hand side
 * \param[in] amg: options of the BoomerAMG preconditioner
 */
static int solve(mfem::ParFiniteElementSpace& fespace,
                 mfem::HypreParMatrix& A,
                 const mfem::Vector& B,
                 const mfem_mgis::Parameters& amg) {
  using namespace mfem_mgis;
  auto ctx = Context{};
  auto& factory = LinearSolverFactory<true>::getFactory();
  auto h = factory.generate(
      ctx, "HyprePCG", fespace,
      Parameters{{"Tolerance", 1e-10},
                 {"MaximumNumberOfIterations", 1000},
                 {"VerbosityLevel", 0},
                 {"Preconditioner",
                  Parameters{{"Name", "HypreBoomerAMG"}, {"Options", amg}}}});
  if (isInvalid(h)) {
    std::cerr << ctx.getErrorMessage() << std::endl;
    return -1;
  }
  auto n = -1;
  for (int i = 0; i != 2; ++i) {
    auto X = mfem::Vector(B.Size());
    X = 0.0;
    h.linear_solver->SetOperator(A);
    h.linear_solver->Mult(B, X);
    const auto oit = getNumberOfIterationsAtConvergence(*(h.linear_solver));
    if ((!oit.has_value()) || ((i == 1) && (*oit != n))) {
      return -1;
    }
    n = *oit;
  }
  return n;
}  // end of solve

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  auto smesh = mfem::Mesh::MakeCartesian3D(8, 8, 8, mfem::Element::HEXAHEDRON);
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  auto fec = mfem::H1_FECollection(1, 3);
  // unknowns ordered by nodes, as in mfem-mgis
  auto fespace = mfem::ParFiniteElementSpace(&pmesh, &fec, 3);
  // linear elasticity, clamped on the first boundary
  auto lambda = mfem::ConstantCoefficient(1.5);
  auto mu = mfem::ConstantCoefficient(1);
  auto a = mfem::ParBilinearForm(&fespace);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
  a.Assemble();
  auto f = mfem::Vector(3);
  f = 0.0;
  f(2) = -1;
  auto fc = mfem::VectorConstantCoefficient(f);
  auto b = mfem::ParLinearForm(&fespace);
  b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(fc));
  b.Assemble();
  auto x = mfem::ParGridFunction(&fespace);
  x = 0.0;
  auto boundaries = mfem::Array<int>(pmesh.bdr_attributes.Max());
  boundaries = 0;
  boundaries[0] = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
  //
  using mfem_mgis::Parameters;
  const auto n_scalar = solve(
      fespace, A, B, Parameters{{"VerbosityLevel", 0}, {"Strategy", "None"}});
  const auto n_default =
      solve(fespace, A, B, Parameters{{"VerbosityLevel", 0}});
  const auto n_system = solve(
      fespace, A, B, Parameters{{"VerbosityLevel", 0}, {"Strategy", "System"}});
  if (fespace.GetMyRank() == 0) {
    std::cout << "number of iterations: " << n_scalar << " (scalar AMG), "
              << n_default << " (default), " << n_system << " (system AMG)\n";
  }
  // a system AMG must take clearly fewer iterations than a scalar AMG
  const auto success = (n_scalar > 0) && (n_default > 0) && (n_system > 0) &&
                       (3 * n_default < 2 * n_scalar) &&
                       (3 * n_system < 2 * n_scalar);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
