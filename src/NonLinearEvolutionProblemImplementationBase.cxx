/*!
 * \file   src/NonLinearEvolutionProblemImplementationBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  const char* const NonLinearEvolutionProblemImplementationBase::
      UseMultiMaterialNonLinearIntegrator =
          "UseMultiMaterialNonLinearIntegrator";

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
      mgis::raise(msg);
    }
  }  // end of checkMultiMaterialSupportEnabled

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getMaterialIdentifiers() const {
    checkMultiMaterialSupportEnabled("getMaterialIdentifiers",
                                     this->mgis_integrator);
    return this->mgis_integrator->getMaterialIdentifiers();
  }  // end of getMaterialIdentifiers

  void NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    checkMultiMaterialSupportEnabled("addBehaviourIntegrator",
                                     this->mgis_integrator);
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
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

  void NonLinearEvolutionProblemImplementationBase::setup(const real t,
                                                          const real dt) {
    if (this->initialization_phase) {
      if (!this->dirichlet_boundary_conditions.empty()) {
        auto fixed_dirichlet_dofs = std::vector<mfem_mgis::size_type>{};
        for (const auto& bc : this->dirichlet_boundary_conditions) {
          auto dofs = bc->getHandledDegreesOfFreedom();
          fixed_dirichlet_dofs.insert(fixed_dirichlet_dofs.end(), dofs.begin(),
                                      dofs.end());
        }
        this->markDegreesOfFreedomHandledByDirichletBoundaryConditions(
            fixed_dirichlet_dofs);
      }
    }
    this->initialization_phase = false;
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      bc->updateImposedValues(this->u1, t + dt);
    }
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setup(t, dt);
    }
  }  // end of NonLinearEvolutionProblemImplementationBase::setup

  void NonLinearEvolutionProblemImplementationBase::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of
     // NonLinearEvolutionProblemImplementationBase::addBoundaryCondition

  NonLinearEvolutionProblemImplementationBase::
      ~NonLinearEvolutionProblemImplementationBase() = default;

}  // end of namespace mfem_mgis
