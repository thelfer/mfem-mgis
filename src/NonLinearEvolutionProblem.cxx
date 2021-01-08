/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/ResidualOperator.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : mfem::NonlinearForm(&(fed->getFiniteElementSpace())),
        mgis_integrator(new MultiMaterialNonLinearIntegrator(fed, h)),
        fe_discretization(fed),
        hypothesis(h),
        u0(fed->getFiniteElementSpace().GetTrueVSize()),
        u1(fed->getFiniteElementSpace().GetTrueVSize()) {
    this->residual = std::make_unique<ResidualOperator>(*this);
    this->u0 = real{0};
    this->u1 = real{0};
    this->solver.SetOperator(*(this->residual));
    this->solver.iterative_mode = true;
    if (this->fe_discretization->getMesh().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
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
    return this->fe_discretization->getFiniteElementSpace();
  }  // end of NonLinearEvolutionProblem::getFiniteElementSpace

  mfem::FiniteElementSpace& NonLinearEvolutionProblem::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace();
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
    this->mgis_integrator->setTimeIncrement(dt);
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  NonLinearEvolutionProblem::~NonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis
