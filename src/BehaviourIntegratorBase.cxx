/*!
 * \file   src/BehaviourIntegratorBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  BehaviourIntegratorBase::BehaviourIntegratorBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::shared_ptr<const Behaviour> b_ptr)
      : Material(s, b_ptr) {
  }  // end of BehaviourIntegratorBase::BehaviourIntegratorBase

  void BehaviourIntegratorBase::setTimeIncrement(const real dt){
    this->time_increment = dt;
  }  // end of setTimeIncrement

  Material& BehaviourIntegratorBase::getMaterial() {
    return *this;
  }  // end of getMaterial

  const Material& BehaviourIntegratorBase::getMaterial() const {
    return *this;
  }  // end of getMaterial

  void BehaviourIntegratorBase::revert(){
    mgis::behaviour::revert(*this);
  }  // end of revert

  void BehaviourIntegratorBase::update(){
    mgis::behaviour::update(*this);
  }  // end of update

  void BehaviourIntegratorBase::checkHypotheses(const Hypothesis h) const {
    using namespace mgis::behaviour;
    if (this->b.hypothesis != h) {
      const auto h1 = std::string(toString(this->b.hypothesis));
      const auto h2 = std::string(toString(h));
      mgis::raise(
          "BehaviourIntegratorBase::checkHypotheses: "
          "the behaviour hypothesis (" +
          h1 + ") does not match the integrator hypothesis (" + h2 + ")");
    }
  }  // end of BehaviourIntegratorBase::checkHypotheses

  void BehaviourIntegratorBase::integrate(const size_type ip) {
    constexpr const auto it = mgis::behaviour::IntegrationType::
        INTEGRATION_CONSISTENT_TANGENT_OPERATOR;
    const auto r =
        mgis::behaviour::integrate(*this, it, this->time_increment, ip, ip + 1);
    if ((r != 0) && (r != 1)) {
      mgis::raise(
          "BehaviourIntegratorBase::integrate: "
          "behaviour integration failed (" +
          std::to_string(r) + ")");
    }
  }  // end of BehaviourIntegratorBase::integrate

  BehaviourIntegratorBase::~BehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis
