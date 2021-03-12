/*!
 * \file   src/NonLinearEvolutionProblemCommon.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
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

  void NonLinearEvolutionProblemCommon::setup(const real, const real) {
  }  // end of NonLinearEvolutionProblemCommon::setup

  NonLinearEvolutionProblemCommon::~NonLinearEvolutionProblemCommon() = default;

}  // end of namespace mfem_mgis
