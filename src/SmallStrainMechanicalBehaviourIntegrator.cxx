/*!
 * \file   src/SmallStrainMechanicalBehaviourIntegrator.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.hxx"

namespace mfem_mgis {

  static const mfem::IntegrationRule &getIntegrationRule(
      const mfem::FiniteElement &el,
      const mfem::ElementTransformation &Trans,
      const Behaviour &b) {
    const auto order = 2 * Trans.OrderGrad(&el);
    if ((b.hypothesis == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) ||
        (b.hypothesis == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS) ||
        (b.hypothesis == Hypothesis::AXISYMMETRICAL)) {
      const auto& ir = mfem::IntRules.Get(el.GetGeomType(), order + 1);
      return ir;
    }
    const auto& ir = mfem::IntRules.Get(el.GetGeomType(), order);
    return ir;
  } // end of getIntegrationRule

  static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
      const mfem::FiniteElementSpace &fs,
      const size_type m,
      const Behaviour &b) {
    auto selector = [&b](const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans)
        -> const mfem::IntegrationRule & {
      const auto &ir = mfem_mgis::getIntegrationRule(el, Trans, b);
      return ir;
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fs, m, selector);
  }  // end of buildQuadratureSpace

  SmallStrainMechanicalBehaviourIntegratorBase::
      SmallStrainMechanicalBehaviourIntegratorBase(
          const mfem::FiniteElementSpace &fs,
          const size_type m,
          std::shared_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(buildQuadratureSpace(fs, m, *b_ptr), b_ptr) {
  }  // end of SmallStrainMechanicalBehaviourIntegratorBase

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

  const mfem::IntegrationRule &
  SmallStrainMechanicalBehaviourIntegratorBase::getIntegrationRule(
      const mfem::FiniteElement &e,
      const mfem::ElementTransformation &t) const {
    return mfem_mgis::getIntegrationRule(e, t, this->b);
  } // end of getIntegrationRule

  SmallStrainMechanicalBehaviourIntegratorBase::
      ~SmallStrainMechanicalBehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis



