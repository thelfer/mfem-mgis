/*!
 * \file   src/BehaviourIntegrator.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#include "MFEMMGIS/BehaviourIntegrator.hxx"

namespace mfem_mgis{

  BehaviourIntegrator::~BehaviourIntegrator() = default;

} // end of namespace mfem_mgis


  //   const mfem::IntegrationRule &MGISIntegrator::getIntegrationRule(
  //       const mfem::FiniteElement &el, const mfem::ElementTransformation
  //       &Trans) const {
  //     const auto *ir = this->IntRule;
  //     if (ir == nullptr) {
  //       int order = 2 * Trans.OrderGrad(&el);  // correct order?
  //       return mfem::IntRules.Get(el.GetGeomType(), order);
  //     }
  //     return *ir;
  //   }  // end of MGISIntegrator::getIntegrationRule()
