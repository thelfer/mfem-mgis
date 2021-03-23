/*!
 * \file   src/NonLinearEvolutionProblemCommon.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include  <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemCommon.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblemCommon::NonLinearEvolutionProblemCommon(
      std::shared_ptr<FiniteElementDiscretization> fed)
      : fe_discretization(fed),
        u0(getTrueVSize(*fed)),
        u1(getTrueVSize(*fed)) {
    this->u0 = real{0};
    this->u1 = real{0};
  }  // end of NonLinearEvolutionProblemCommon

  FiniteElementDiscretization&
  NonLinearEvolutionProblemCommon::getFiniteElementDiscretization() {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblemCommon::getFiniteElementDiscretizationPointer() {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  mfem::Vector&
  NonLinearEvolutionProblemCommon::getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemCommon::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblemCommon::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblemCommon::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  void NonLinearEvolutionProblemCommon::revert() {
    this->u1 = this->u0;
  }  // end of revert

  void NonLinearEvolutionProblemCommon::update() {
    this->u0 = this->u1;
  }  // end of update

  void NonLinearEvolutionProblemCommon::setTimeIncrement(const real) {
  }  // end of NonLinearEvolutionProblemCommon::setTimeIncrement

  void NonLinearEvolutionProblemCommon::setup(const real t, const real dt) {
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
  }  // end of NonLinearEvolutionProblemCommon::setup

  void NonLinearEvolutionProblemCommon::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->dirichlet_boundary_conditions.push_back(std::move(bc));
  }  // end of NonLinearEvolutionProblemCommon::addBoundaryCondition

  NonLinearEvolutionProblemCommon::~NonLinearEvolutionProblemCommon() = default;

}  // end of namespace mfem_mgis
