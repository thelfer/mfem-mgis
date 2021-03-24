/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    using SequentialImplementation =
        NonLinearEvolutionProblemImplementation<false>;
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      using ParallelImplementation =
          NonLinearEvolutionProblemImplementation<true>;
      this->pimpl = std::make_unique<ParallelImplementation>(fed, h, p);
#else
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      this->pimpl = std::make_unique<SequentialImplementation>(fed, h, p);
    }
  }  // end of NonLinearEvolutionProblem

  FiniteElementDiscretization&
  NonLinearEvolutionProblem::getFiniteElementDiscretization() {
    return this->pimpl->getFiniteElementDiscretization();
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblem::getFiniteElementDiscretizationPointer() {
    return this->pimpl->getFiniteElementDiscretizationPointer();
  }  // end of getFiniteElementDiscretizationPointer

  mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() {
    return this->pimpl->getUnknownsAtBeginningOfTheTimeStep();
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->pimpl->getUnknownsAtBeginningOfTheTimeStep();
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep() {
    return this->pimpl->getUnknownsAtEndOfTheTimeStep();
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->pimpl->getUnknownsAtEndOfTheTimeStep();
  }  // end of getUnknownsAtEndOfTheTimeStep

  NewtonSolver& NonLinearEvolutionProblem::getSolver() {
    return this->pimpl->getSolver();
  }  // end of getSolver

  void NonLinearEvolutionProblem::solve(const real t, const real dt) {
    this->setup(t, dt);
    this->pimpl->solve(t, dt);
  }  // end of solve

  void NonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->pimpl->addBoundaryCondition(std::move(bc));
  }  // end of NonLinearEvolutionProblem::addBoundaryCondition

  void NonLinearEvolutionProblem::addPostProcessing(
      const std::function<void(const real, const real)>& p) {
    this->pimpl->addPostProcessing(p);
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem::executePostProcessings(const real t,
                                                         const real dt) {
    this->pimpl->executePostProcessings(t, dt);
  }

  void NonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
                                                         const size_type l,
                                                         const std::string& m,
                                                         const std::string& b) {
    this->pimpl->addBehaviourIntegrator(n, l, m, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblem::getMaterial(
      const size_type m) const {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblem::getMaterial(const size_type m) {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator& NonLinearEvolutionProblem::getBehaviourIntegrator(
      const size_type m) const {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator& NonLinearEvolutionProblem::getBehaviourIntegrator(
      const size_type m) {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblem::update() {
    this->pimpl->update();
  }  // end of update

  void NonLinearEvolutionProblem::revert() {
    this->pimpl->revert();
  }  // end of revert

  void NonLinearEvolutionProblem::setup(const real, const real) {
  }  // end of setup

  NonLinearEvolutionProblem::~NonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis