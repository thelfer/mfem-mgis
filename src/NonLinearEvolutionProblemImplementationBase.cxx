/*!
 * \file   src/NonLinearEvolutionProblemImplementationBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <iostream>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  const char* const NonLinearEvolutionProblemImplementationBase::
      UseMultiMaterialNonLinearIntegrator =
          "UseMultiMaterialNonLinearIntegrator";

  std::vector<std::string>
  NonLinearEvolutionProblemImplementationBase::getParametersList() {
    return {NonLinearEvolutionProblemImplementationBase::
                UseMultiMaterialNonLinearIntegrator};
  }  // end of getParametersList

  MultiMaterialNonLinearIntegrator* buildMultiMaterialNonLinearIntegrator(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    const auto* const n = NonLinearEvolutionProblemImplementationBase::
        UseMultiMaterialNonLinearIntegrator;
    if (contains(p, n)) {
      if (!get<bool>(p, n)) {
        return nullptr;
      }
    }
    return new MultiMaterialNonLinearIntegrator(fed, h);
  }  // end of buildMultiMaterialNonLinearIntegrator

  NonLinearEvolutionProblemImplementationBase::
      NonLinearEvolutionProblemImplementationBase(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : fe_discretization(fed),
        u0(getTrueVSize(*fed)),
        u1(getTrueVSize(*fed)),
        mgis_integrator(buildMultiMaterialNonLinearIntegrator(fed, h, p)),
        hypothesis(h) {
    this->u0 = real{0};
    this->u1 = real{0};
  }  // end of NonLinearEvolutionProblemImplementationBase

  FiniteElementDiscretization& NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretization() {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  const FiniteElementDiscretization&
  NonLinearEvolutionProblemImplementationBase::getFiniteElementDiscretization()
      const {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretizationPointer() {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  void NonLinearEvolutionProblemImplementationBase::revert() {
    this->u1 = this->u0;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->revert();
    }
  }  // end of revert

  void NonLinearEvolutionProblemImplementationBase::update() {
    this->u0 = this->u1;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->update();
    }
  }  // end of update

  static void checkMultiMaterialSupportEnabled(
      const char* const n, const MultiMaterialNonLinearIntegrator* const p) {
    if (p == nullptr) {
      std::string msg("NonLinearEvolutionProblemImplementationBase::");
      msg += n;
      msg += ": multi material support has been disabled";
      raise(msg);
    }
  }  // end of checkMultiMaterialSupportEnabled

  void NonLinearEvolutionProblemImplementationBase::setMaterialsNames(
      const std::map<size_type, std::string>& ids) {
    this->getFiniteElementDiscretization().setMaterialsNames(ids);
  }

  void NonLinearEvolutionProblemImplementationBase::setBoundariesNames(
      const std::map<size_type, std::string>& ids){
    this->getFiniteElementDiscretization().setBoundariesNames(ids);
  }

  size_type NonLinearEvolutionProblemImplementationBase::getMaterialIdentifier(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getMaterialIdentifier(p);
  }  // end of getMaterialIdentifier

  size_type NonLinearEvolutionProblemImplementationBase::getBoundaryIdentifier(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getBoundaryIdentifier(p);
  }  // end of getBoundaryIdentifier

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getMaterialsIdentifiers(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getMaterialsIdentifiers(p);
  }  // end of getMaterialsIdentifiers

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getAssignedMaterialsIdentifiers() const {
    checkMultiMaterialSupportEnabled("getAssignedMaterialsIdentifiers",
                                     this->mgis_integrator);
    return this->mgis_integrator->getAssignedMaterialsIdentifiers();
  }  // end of getAssignedMaterialsIdentifiers

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getBoundariesIdentifiers(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getBoundariesIdentifiers(p);
  }  // end of getBoundariesIdentifiers

  void NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      const std::string& n,
      const Parameter& m,
      const std::string& l,
      const std::string& b) {
    checkMultiMaterialSupportEnabled("addBehaviourIntegrator",
                                     this->mgis_integrator);
    for (const auto& id : this->getMaterialsIdentifiers(m)) {
      this->mgis_integrator->addBehaviourIntegrator(n, id, l, b);
    }
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const size_type m) const {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const size_type m) {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) const {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblemImplementationBase::setTimeIncrement(
      const real dt) {
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setTimeIncrement(dt);
    }
  }  // end of setTimeIncrement

  void NonLinearEvolutionProblemImplementationBase::setMacroscopicGradients(
      const std::vector<real>& g) {
    checkMultiMaterialSupportEnabled("setMacroscopicGradients",
                                     this->mgis_integrator);
    this->mgis_integrator->setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  void NonLinearEvolutionProblemImplementationBase::setSolverParameters(
      const Parameters& params) {
#ifdef MFEM_USE_PETSC
    if (usePETSc()) {
      mfem_mgis::setSolverParameters(*(this->petsc_solver), params);
    } else {
      mfem_mgis::setSolverParameters(*(this->solver), params);
    }
#else  /* MFEM_USE_PETSC */
    mfem_mgis::setSolverParameters(*(this->solver), params);
#endif /* MFEM_USE_PETSC */
  }  // end of setSolverParameters

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getEssentialDegreesOfFreedom()
      const {
    auto ddofs = std::vector<mfem_mgis::size_type>{};
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      auto dofs = bc->getHandledDegreesOfFreedom();
      ddofs.insert(ddofs.end(), dofs.begin(), dofs.end());
    }
    return ddofs;
  }  // end of getEssentialDegreesOfFreedom

  void NonLinearEvolutionProblemImplementationBase::setup(const real t,
                                                          const real dt) {
    if (this->initialization_phase) {
      if (!this->dirichlet_boundary_conditions.empty()) {
        this->markDegreesOfFreedomHandledByDirichletBoundaryConditions(
            this->getEssentialDegreesOfFreedom());
      }
    }
    this->initialization_phase = false;
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      bc->updateImposedValues(this->u1, t + dt);
    }
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setup(t, dt);
    }
  }  // end of setup

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s) {
    if (usePETSc()) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::updateLinearSolver: "
          "call to this method is meaningless if PETSc is used");
    }
    this->linear_solver_preconditioner.reset();
    this->linear_solver = std::move(s);
    this->solver->setLinearSolver(*(this->linear_solver));
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s,
      std::unique_ptr<LinearSolverPreconditioner> p) {
    if (usePETSc()) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::updateLinearSolver: "
          "call to this method is meaningless if PETSc is used");
    }
    if (p != nullptr) {
      auto* const isolver = dynamic_cast<IterativeSolver*>(s.get());
      if (isolver != nullptr) {
	isolver->SetPreconditioner(*p);
      }
      this->updateLinearSolver(std::move(s));
      this->linear_solver_preconditioner = std::move(p);
    } else {
      this->updateLinearSolver(std::move(s));
    }
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      LinearSolverHandler s) {
    this->updateLinearSolver(std::move(s.linear_solver),
                             std::move(s.preconditioner));
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblemImplementationBase::solve(const real t,
                                                          const real dt) {
    this->setTimeIncrement(dt);
    this->setup(t, dt);
    //    this->computePrediction(t, dt);
    auto success = false;
    if (usePETSc()) {
#ifdef MFEM_USE_PETSC
      mfem::Vector zero;
      // PETSc solver somehow "requires" zero as a first argument of Mult.
      // If it is not the case, one should take care of memory management.
      this->petsc_solver->Mult(zero, this->u1);
      success = this->petsc_solver->GetConverged();
#else  /* MFEM_USE_PETSC */
      MFEM_VERIFY(0, "Support for PETSc is deactivated");
#endif /* MFEM_USE_PETSC */
    } else {
      this->solver->Mult(this->u0, this->u1);
      success = this->solver->GetConverged();
    }
    return success;
  }  // end of solve

  void NonLinearEvolutionProblemImplementationBase::computePrediction(const real /*t*/, const real /*dt*/) {
    // if (mgis_integrator == nullptr) {
    //   return;
    // }
    // constexpr auto it = IntegrationType::PREDICTION_ELASTIC_OPERATOR;
    // if (!this->integrate(this->u0, it)) {
    //   return;
    // }
    // mfem::Vector du(this->u1.Size());
    // mfem::Vector r(this->u1.Size());
    // r = 0.;
    // du = 0.;
    // //
    // // unmarked Dirichlet boundary conditions
    // const auto ddofs = this->getEssentialDegreesOfFreedom();
    // this->markDegreesOfFreedomHandledByDirichletBoundaryConditions({});
    // auto& K = this->solver->getJacobian(this->u0);
    // for (const auto& bc : this->dirichlet_boundary_conditions) {
    //   bc->setImposedValuesIncrements(du, t, t + dt);
    // }

    // // solve the tangent problem
    // this->markDegreesOfFreedomHandledByDirichletBoundaryConditions(ddofs);
    // for (const auto& bc : this->dirichlet_boundary_conditions) {
    //   bc->setImposedValuesIncrements(r, t, t + dt);
    // }
    // if (this->solver->computeNewtonCorrection(du, r, this->u0)) {
    //   this->u1 = this->u0;
    //   this->u1 += du;
    // }
  }  // end of computePrediction

  NonLinearEvolutionProblemImplementationBase::
      ~NonLinearEvolutionProblemImplementationBase() = default;

}  // end of namespace mfem_mgis
