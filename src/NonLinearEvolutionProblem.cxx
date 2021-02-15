/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  NonLinearEvolutionProblem<true>::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : NonLinearEvolutionProblemBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h) {
    if (this->fe_discretization->getMesh<true>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->AddDomainIntegrator(this->parallel_mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem<true>::setup() {
    MultiMaterialEvolutionProblemBase::setup();
  }  // end of setup

  void NonLinearEvolutionProblem<true>::revert() {
    NonLinearEvolutionProblemBase<true>::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblem<true>::update() {
    NonLinearEvolutionProblemBase<false>::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblem<true>::setTimeIncrement(const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

  NonLinearEvolutionProblem<true>::~NonLinearEvolutionProblem() = default;

#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblem<false>::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : NonLinearEvolutionProblemBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h) {
    if (this->fe_discretization->getMesh<false>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->AddDomainIntegrator(this->sequential_mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem<false>::setup() {
    MultiMaterialEvolutionProblemBase::setup();
  }  // end of setup

  void NonLinearEvolutionProblem<false>::revert() {
    NonLinearEvolutionProblemBase<false>::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblem<false>::update() {
    NonLinearEvolutionProblemBase<false>::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblem<false>::setTimeIncrement(const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

  NonLinearEvolutionProblem<false>::~NonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis
