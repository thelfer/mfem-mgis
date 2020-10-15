/*!
 * \file   src/SmallStrainMechanicalBehaviourIntegrator.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
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
      return mfem::IntRules.Get(el.GetGeomType(), order + 1);
    }
    return mfem::IntRules.Get(el.GetGeomType(), order);
  } // end of getIntegrationRule

  static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
      const mfem::FiniteElementSpace &fs,
      const size_type m,
      const Behaviour &b) {
    auto selector = [&b](const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans) {
      return getIntegrationRule(el, Trans, b);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fs, m, selector);
  }  // end of buildQuadratureSpace

  SmallStrainMechanicalBehaviourIntegratorBase::
      SmallStrainMechanicalBehaviourIntegratorBase(
          const mfem::FiniteElementSpace &fs,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(buildQuadratureSpace(fs, m, *b_ptr),
                                std::move(b_ptr)) {
  }  // end of SmallStrainMechanicalBehaviourIntegrator

  SmallStrainMechanicalBehaviourIntegratorBase::
      ~SmallStrainMechanicalBehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis



