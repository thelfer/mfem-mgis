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
        hypothesis(h),
        u0(fs->GetTrueVSize()),
        u1(fs->GetTrueVSize()) {
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::FiniteElementSpace&
  NonLinearEvolutionProblem::getFiniteElementSpace() const {
    return *(this->fe_space);
  }  // end of NonLinearEvolutionProblem::getFiniteElementSpace

  mfem::NewtonSolver& NonLinearEvolutionProblem::getSolver() {
    return this->solver;
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem::revert() {
    this->u1 = this->u0;
    this->mgis_integrator->revert();
  }  // end of revert

  void NonLinearEvolutionProblem::update() {
    this->u0 = this->u1;
    this->mgis_integrator->update();
  }  // end of update

  void NonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
                                                         const size_type m,
                                                         const std::string& l,
                                                         const std::string& b) {
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblem::getMaterial(
      const size_type m) const {
    return this->mgis_integrator->getMaterial(m);
  }  // end of NonLinearEvolutionProblem::getMaterial

  Material& NonLinearEvolutionProblem::getMaterial(const size_type m) {
    return this->mgis_integrator->getMaterial(m);
  }  // end of NonLinearEvolutionProblem::getMaterial

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
