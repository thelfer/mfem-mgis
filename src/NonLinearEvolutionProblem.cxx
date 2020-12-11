/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/ResidualOperator.hxx"
#include "MFEMMGIS/MGISIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<mfem::FiniteElementSpace> fs, const Hypothesis h)
      : mfem::NonlinearForm(fs.get()),
        mgis_integrator(new MGISIntegrator(fs, h)),
        fe_space(fs),
        hypothesis(h) {
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  mfem::NewtonSolver& NonLinearEvolutionProblem::getSolver() {
    return this->solver;
    }  // end of NonLinearEvolutionProblem

    void NonLinearEvolutionProblem::solve(const real dt) {
      mfem::Vector zero;
      this->u1 = this->u0;
      this->mgis_integrator->setTimeIncrement(dt);
      ResidualOperator r(*this);
      this->solver.SetOperator(r);
      this->solver.Mult(zero, this->u1);
      if (!this->solver.GetConverged()) {
        mgis::raise("Newton solver did not converge");
      }
    }  // end of solve

    NonLinearEvolutionProblem::~NonLinearEvolutionProblem() = default;

  }  // end of namespace mfem_mgis
