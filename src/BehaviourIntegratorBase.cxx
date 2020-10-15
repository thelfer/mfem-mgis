/*!
 * \file   src/BehaviourIntegratorBase.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis{

  BehaviourIntegratorBase::BehaviourIntegratorBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::unique_ptr<const Behaviour> b_ptr)
      : Material(s, std::move(b_ptr))
  {}  // end of BehaviourIntegratorBase::BehaviourIntegratorBase

  void BehaviourIntegratorBase::checkHypotheses(const Hypothesis h) const {
    using namespace mgis::behaviour;
    if (this->b.hypothesis != h) {
      const auto h1 = std::string(toString(this->b.hypothesis));
      const auto h2 = std::string(toString(h));
      mgis::raise(
          "SmallStrainMechanicalBehaviourIntegrator::"
          "SmallStrainMechanicalBehaviourIntegrator: "
          "the behaviour hypothesis (" +
          h1 + ") does not match the integrator hypothesis (" + h2 + ")");
    }
  }  // end of BehaviourIntegratorBase::checkHypotheses

  BehaviourIntegratorBase::~BehaviourIntegratorBase() = default;

} // end of namespace mfem_mgis
