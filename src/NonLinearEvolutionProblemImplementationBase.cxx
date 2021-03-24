/*!
 * \file   src/NonLinearEvolutionProblemImplementationBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblemImplementationBase::
      NonLinearEvolutionProblemImplementationBase(
          std::shared_ptr<FiniteElementDiscretization> fed)
      : fe_discretization(fed), u0(getTrueVSize(*fed)), u1(getTrueVSize(*fed)) {
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

  mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  void NonLinearEvolutionProblemImplementationBase::revert() {
    this->u1 = this->u0;
  }  // end of revert

  void NonLinearEvolutionProblemImplementationBase::update() {
    this->u0 = this->u1;
  }  // end of update

  void NonLinearEvolutionProblemImplementationBase::setTimeIncrement(
      const real) {
  }  // end of NonLinearEvolutionProblemImplementationBase::setTimeIncrement

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
  }  // end of NonLinearEvolutionProblemImplementationBase::setup

  void NonLinearEvolutionProblemImplementationBase::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of
     // NonLinearEvolutionProblemImplementationBase::addBoundaryCondition

  NonLinearEvolutionProblemImplementationBase::
      ~NonLinearEvolutionProblemImplementationBase() = default;

}  // end of namespace mfem_mgis
