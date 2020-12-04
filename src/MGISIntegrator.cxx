/*!
 * \file   src/MGISIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/MGISIntegrator.hxx"

namespace mfem_mgis {

  MGISIntegrator::MGISIntegrator(
      std::shared_ptr<const mfem::FiniteElementSpace> fs, const Hypothesis h)
      : fe_space(std::move(fs)),
        hypothesis(h) {}  // end of MGISIntegrator::MGISIntegrator

  void MGISIntegrator::AssembleElementVector(const mfem::FiniteElement & e,
                                             mfem::ElementTransformation & tr,
                                             const mfem::Vector & U,
                                             mfem::Vector & F) {
    const auto m = tr.Attribute;
    const auto pbi = this->behaviour_integrators.find(m);
    if (pbi == this->behaviour_integrators.end()) {
      mgis::raise(
          "MGISIntegrator::AssembleElementVector: "
          "no behaviour integrator associated with material '" +
          std::to_string(m) + "'");
    }
    pbi->second->computeInnerForces(F, e, tr, U);
  }  // end of MGISIntegrator::AssembleElementVector

  void MGISIntegrator::AssembleElementGrad(const mfem::FiniteElement& e,
                                           mfem::ElementTransformation& tr,
                                           const mfem::Vector& U,
                                           mfem::DenseMatrix& K) {
    const auto m = tr.Attribute;
    const auto pbi = this->behaviour_integrators.find(m);
    if (pbi == this->behaviour_integrators.end()) {
      mgis::raise(
          "MGISIntegrator::AssembleElementGrad: "
          "no behaviour integrator associated with material '" +
          std::to_string(m) + "'");
    }
    pbi->second->computeStiffnessMatrix(K, e, tr, U);
  }  // end of MGISIntegrator::AssembleElementGrad

  void MGISIntegrator::addBehaviourIntegrator(const std::string& n,
                                              const size_type m,
                                              const std::string& l,
                                              const std::string& b) {
    if (this->behaviour_integrators.count(m) != 0) {
      mgis::raise(
          "MGISIntegrator::addBehaviourIntegrator: "
          "integrator already defined for material '" +
          std::to_string(m) + "'");
    }
    const auto& f = BehaviourIntegratorFactory::get(this->hypothesis);
    this->behaviour_integrators[m] = f.generate(
        n, *(this->fe_space), m, mfem_mgis::load(l, b, this->hypothesis));
  }  // end of MGISIntegrator::addBehaviourIntegrator

  void MGISIntegrator::setTimeIncrement(const real dt) {
    this->time_increment = dt;
  } // end of MGISIntegrator::setTimeIncrement

  MGISIntegrator::~MGISIntegrator() = default;

}  // end of namespace mfem_mgis
