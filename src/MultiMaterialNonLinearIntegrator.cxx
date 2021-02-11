/*!
 * \file   src/MultiMaterialNonLinearIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple test
   * \param[in] i: pointer to behaviour integrator
   * \param[in] n: name of the calling method
   * \param[in] m: material id
   */
  static void checkIfBehaviourIntegratorIsDefined(
      const BehaviourIntegrator* const i,
      const char* const n,
      const size_type m) {
    if (i == nullptr) {
      mgis::raise("MultiMaterialNonLinearIntegrator::" + std::string(n) +
                  ": no behaviour integrator associated with material '" +
                  std::to_string(m) + "'");
    }
  }  // end if checkIfBehaviourIntegratorIsDefined

  MultiMaterialNonLinearIntegrator::MultiMaterialNonLinearIntegrator(
      std::shared_ptr<const FiniteElementDiscretization> fed, const Hypothesis h)
      : fe_discretization(fed),
        hypothesis(h) {
    const auto& mesh = this->fe_discretization->getMesh();
    // shifting by one allows to directly use the material id to get
    // the behaviour integrators in all the other methods.
    // However, behaviour_integrators[0] will always be null.
    this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
  }  // end of MultiMaterialNonLinearIntegrator

  void MultiMaterialNonLinearIntegrator::AssembleElementVector(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::Vector& F) {
    const auto m = tr.Attribute;
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementVector", m);
    bi->computeInnerForces(F, e, tr, U);
  }  // end of AssembleElementVector

  void MultiMaterialNonLinearIntegrator::AssembleElementGrad(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::DenseMatrix& K) {
    const auto m = tr.Attribute;
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementGrad", m);
    bi->computeStiffnessMatrix(K, e, tr, U);
  }  // end of AssembleElementGrad

  void MultiMaterialNonLinearIntegrator::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    if (this->behaviour_integrators[m] != nullptr) {
      mgis::raise(
          "MultiMaterialNonLinearIntegrator::addBehaviourIntegrator: "
          "integrator already defined for material '" +
          std::to_string(m) + "'");
    }
    const auto& f = BehaviourIntegratorFactory::get(this->hypothesis);
    this->behaviour_integrators[m] =
        f.generate(n, *(this->fe_discretization), m,
                   mfem_mgis::load(l, b, this->hypothesis));
  }  // end of addBehaviourIntegrator

  const Material& MultiMaterialNonLinearIntegrator::getMaterial(
      const size_type m) const {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
    return bi->getMaterial();
  }  // end of getMaterial

  Material& MultiMaterialNonLinearIntegrator::getMaterial(const size_type m) {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
    return bi->getMaterial();
  }  // end of getMaterial

  void MultiMaterialNonLinearIntegrator::setTimeIncrement(const real dt) {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->setTimeIncrement(dt);
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setup() {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->setup();
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::revert() {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->revert();
      }
    }
  }  // end of revert

  void MultiMaterialNonLinearIntegrator::update() {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->update();
      }
    }
  }  // end of update

  MultiMaterialNonLinearIntegrator::~MultiMaterialNonLinearIntegrator() =
      default;

}  // end of namespace mfem_mgis
