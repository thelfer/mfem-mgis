/*!
 * \file   src/MGISIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/MGISIntegrator.hxx"

namespace mfem_mgis {

  MGISIntegrator::MGISIntegrator(
      std::shared_ptr<const mfem::FiniteElementSpace> fs)
      : fe_space(std::move(fs)) {}  // end of MGISIntegrator::MGISIntegrator

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
    pbi->second->computeInnerForces(e, tr, U, F);
  }  // end of MGISIntegrator::AssembleElementVector

  void MGISIntegrator::AssembleElementGrad(const mfem::FiniteElement& e,
                                           mfem::ElementTransformation& tr,
                                           const mfem::Vector& U,
                                           mfem::DenseMatrix& K) {
    const auto m = tr.Attribute;
    const auto pbi = this->behaviour_integrators.find(m);
    if (pbi == this->behaviour_integrators.end()) {
      mgis::raise(
          "MGISIntegrator::AssembleElementVector: "
          "no behaviour integrator associated with material '" +
          std::to_string(m) + "'");
    }
    pbi->second->computeStiffnessMatrix(e, tr, U, K);
  }  // end of MGISIntegrator::AssembleElementGrad

  MGISIntegrator::~MGISIntegrator() = default;

}  // end of namespace mfem_mgis
