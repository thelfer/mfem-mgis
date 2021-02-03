/*!
 * \file   src/NonLinearEvolutionProblemBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemBase.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase(
      std::shared_ptr<FiniteElementDiscretization> fed)
      : mfem::NonlinearForm(&(fed->getFiniteElementSpace())),
        fe_discretization(fed),
        u0(fed->getFiniteElementSpace().GetTrueVSize()),
        u1(fed->getFiniteElementSpace().GetTrueVSize()) {
    this->u0 = real{0};
    this->u1 = real{0};
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
  }  // end of NonLinearEvolutionProblemBase

  mfem::Vector&
  NonLinearEvolutionProblemBase::getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemBase::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblemBase::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblemBase::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::FiniteElementSpace&
  NonLinearEvolutionProblemBase::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace();
  }  // end of NonLinearEvolutionProblemBase::getFiniteElementSpace

  mfem::FiniteElementSpace& NonLinearEvolutionProblemBase::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace();
  }  // end of NonLinearEvolutionProblemBase::getFiniteElementSpace

  mfem::NewtonSolver& NonLinearEvolutionProblemBase::getSolver() {
    return this->solver;
  }  // end of NonLinearEvolutionProblemBase

  void NonLinearEvolutionProblemBase::revert() {
    this->u1 = this->u0;
  }  // end of revert

  void NonLinearEvolutionProblemBase::update() {
    this->u0 = this->u1;
  }  // end of update

  void NonLinearEvolutionProblemBase::setTimeIncrement(const real) {
  }  // end of NonLinearEvolutionProblemBase::setTimeIncrement

  void NonLinearEvolutionProblemBase::solve(const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  NonLinearEvolutionProblemBase::~NonLinearEvolutionProblemBase() = default;

}  // end of namespace mfem_mgis
