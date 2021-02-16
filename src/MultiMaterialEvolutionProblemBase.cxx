/*!
 * \file   MultiMaterialEvolutionProblemBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/MultiMaterialEvolutionProblemBase.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static MultiMaterialNonLinearIntegrator<parallel>*
  buildMultiMaterialNonLinearIntegrator(
      const std::shared_ptr<FiniteElementDiscretization>& fed,
      const Hypothesis h) {
    if constexpr (parallel) {
      if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        return new MultiMaterialNonLinearIntegrator<true>(fed, h);
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      } else {
        return nullptr;
      }
    } else {
      if (fed->describesAParallelComputation()) {
        return nullptr;
      } else {
        return new MultiMaterialNonLinearIntegrator<false>(fed, h);
      }
    }
  }  // end of buildMultiMaterialNonLinearIntegrator

  MultiMaterialEvolutionProblemBase::MultiMaterialEvolutionProblemBase(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      :
#ifdef MFEM_USE_MPI
        parallel_mgis_integrator(
            buildMultiMaterialNonLinearIntegrator<true>(fed, h)),
#endif /* MFEM_USE_MPI */
        sequential_mgis_integrator(
            buildMultiMaterialNonLinearIntegrator<false>(fed, h)),
        hypothesis(h) {
  }  // end of MultiMaterialEvolutionProblemBase

  MultiMaterialNonLinearIntegratorBase&
  MultiMaterialEvolutionProblemBase::getMultiMaterialNonLinearIntegratorBase() {
#ifdef MFEM_USE_MPI
    if (this->sequential_mgis_integrator != nullptr) {
      return *(this->sequential_mgis_integrator);
    }
    return *(this->parallel_mgis_integrator);
#else  /* MFEM_USE_MPI */
    return *(this->sequential_mgis_integrator);
#endif /* MFEM_USE_MPI */
  }    // end of getMultiMaterialNonLinearIntegratorBase

  const MultiMaterialNonLinearIntegratorBase&
  MultiMaterialEvolutionProblemBase::getMultiMaterialNonLinearIntegratorBase()
      const {
#ifdef MFEM_USE_MPI
    if (this->sequential_mgis_integrator != nullptr) {
      return *(this->sequential_mgis_integrator);
    }
    return *(this->parallel_mgis_integrator);
#else  /* MFEM_USE_MPI */
    return *(this->sequential_mgis_integrator);
#endif /* MFEM_USE_MPI */
  }    // end of getMultiMaterialNonLinearIntegratorBase

  void MultiMaterialEvolutionProblemBase::setup() {
    this->getMultiMaterialNonLinearIntegratorBase().setup();
  }  // end of setup

  void MultiMaterialEvolutionProblemBase::revert() {
    this->getMultiMaterialNonLinearIntegratorBase().revert();
  }  // end of revert

  void MultiMaterialEvolutionProblemBase::update() {
    this->getMultiMaterialNonLinearIntegratorBase().update();
  }  // end of update

  void MultiMaterialEvolutionProblemBase::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    this->getMultiMaterialNonLinearIntegratorBase().addBehaviourIntegrator(
        n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& MultiMaterialEvolutionProblemBase::getMaterial(
      const size_type m) const {
    return this->getMultiMaterialNonLinearIntegratorBase().getMaterial(m);
  }  // end of MultiMaterialEvolutionProblemBase::getMaterial

  Material& MultiMaterialEvolutionProblemBase::getMaterial(const size_type m) {
    return this->getMultiMaterialNonLinearIntegratorBase().getMaterial(m);
  }  // end of MultiMaterialEvolutionProblemBase::getMaterial

  void MultiMaterialEvolutionProblemBase::setTimeIncrement(const real dt) {
    this->getMultiMaterialNonLinearIntegratorBase().setTimeIncrement(dt);
  }  // end of setTimeIncrement

  MultiMaterialEvolutionProblemBase::~MultiMaterialEvolutionProblemBase() =
      default;

}  // end of namespace mfem_mgis
