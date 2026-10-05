/*!
 * \file   HypreBoomerAMGTest.cxx
 * \brief  Test of the system and elasticity BoomerAMG preconditioners built
 * by the linear solver factory
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

//! \brief numbers of iterations obtained with the strategies of BoomerAMG
struct NumbersOfIterations {
  //! \brief strategy `None`
  int scalar = -1;
  //! \brief no strategy given
  int by_default = -1;
  //! \brief strategy `System`
  int system = -1;
  //! \brief strategy `Elasticity`
  int elasticity = -1;
};

/*!
 * \return the numbers of iterations obtained with the strategies of BoomerAMG
 * on a cantilever beam in linear elasticity
 * \param[in] pmesh: mesh
 * \param[in] ordering: ordering of the unknowns
 */
static NumbersOfIterations getNumbersOfIterations(
    mfem::ParMesh& pmesh, const mfem::Ordering::Type ordering) {
  auto fec = mfem::H1_FECollection(1, 3);
  auto fespace = mfem::ParFiniteElementSpace(&pmesh, &fec, 3, ordering);
  // linear elasticity, clamped on the face x = 0
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
  boundaries[4] = 1;
  auto ess_tdofs = mfem::Array<int>{};
  fespace.GetEssentialTrueDofs(boundaries, ess_tdofs);
  auto A = mfem::HypreParMatrix{};
  auto B = mfem::Vector{};
  auto X = mfem::Vector{};
  a.FormLinearSystem(ess_tdofs, x, b, A, X, B);
  //
  using mfem_mgis::Parameters;
  auto r = NumbersOfIterations{};
  r.scalar = solve(fespace, A, B,
                   Parameters{{"VerbosityLevel", 0}, {"Strategy", "None"}});
  r.by_default = solve(fespace, A, B, Parameters{{"VerbosityLevel", 0}});
  r.system = solve(fespace, A, B,
                   Parameters{{"VerbosityLevel", 0}, {"Strategy", "System"}});
  r.elasticity =
      solve(fespace, A, B,
            Parameters{{"VerbosityLevel", 0}, {"Strategy", "Elasticity"}});
  if (fespace.GetMyRank() == 0) {
    std::cout << "number of iterations ("
              << (ordering == mfem::Ordering::byNODES ? "byNODES" : "byVDIM")
              << "): " << r.scalar << " (scalar AMG), " << r.by_default
              << " (default), " << r.system << " (system AMG), " << r.elasticity
              << " (elasticity AMG)\n";
  }
  return r;
}  // end of getNumbersOfIterations

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  // slender beam, for which the rigid body modes used by the `Elasticity`
  // strategy matter
  auto smesh =
      mfem::Mesh::MakeCartesian3D(32, 4, 4, mfem::Element::HEXAHEDRON, 8, 1, 1);
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, smesh);
  // a system AMG must take clearly fewer iterations than a scalar AMG
  auto check = [](const NumbersOfIterations& r) {
    return (r.scalar > 0) && (r.by_default > 0) && (r.system > 0) &&
           (r.elasticity > 0) && (3 * r.by_default < 2 * r.scalar) &&
           (3 * r.system < 2 * r.scalar) && (3 * r.elasticity < 2 * r.scalar);
  };
  // unknowns ordered by nodes, the default in mfem-mgis: the `Elasticity`
  // strategy falls back to the `System` strategy
  const auto r1 = getNumbersOfIterations(pmesh, mfem::Ordering::byNODES);
  // unknowns ordered by vector dimension: the `Elasticity` strategy is really
  // used, its number of iterations differs from the one of the `System`
  // strategy
  const auto r2 = getNumbersOfIterations(pmesh, mfem::Ordering::byVDIM);
  const auto success = check(r1) && (r1.elasticity == r1.system) && check(r2) &&
                       (r2.elasticity != r2.system);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
