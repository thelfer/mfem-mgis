/*!
 * \file   src/BehaviourIntegratorBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  BehaviourIntegratorBase::BehaviourIntegratorBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::unique_ptr<const Behaviour> b_ptr)
      : Material(s, std::move(b_ptr)) {
    // The following arrays are storing material properties and external
    // state variables. They can be allocated a single time
    this->wks.mps.resize(getArraySize(this->b.mps, this->b.hypothesis));
    this->wks.esvs0.resize(getArraySize(this->b.esvs, this->b.hypothesis));
    this->wks.esvs1.resize(getArraySize(this->b.esvs, this->b.hypothesis));
  }  // end of BehaviourIntegratorBase

  void BehaviourIntegratorBase::throwInvalidBehaviourType(
      const char* const mn, const char* const m) const {
    auto msg = std::string(mn) + ": invalid behaviour type";
    if (m != nullptr) {
      msg += "(" + std::string(m) + ')';
    }
    mgis::raise(msg);
  }  // end of throwInvalidBehaviourType

  void BehaviourIntegratorBase::throwInvalidBehaviourKinematic(
      const char* const mn, const char* const m) const {
    auto msg = std::string(mn) + ": invalid behaviour type";
    if (m != nullptr) {
      msg += "(" + std::string(m) + ')';
    }
    mgis::raise(msg);
  }  // end of throwInvalidBehaviourType

  void BehaviourIntegratorBase::setTimeIncrement(const real dt){
    this->time_increment = dt;
  }  // end of setTimeIncrement

  Material& BehaviourIntegratorBase::getMaterial() {
    return *this;
  }  // end of getMaterial

  const Material& BehaviourIntegratorBase::getMaterial() const {
    return *this;
  }  // end of getMaterial

  void BehaviourIntegratorBase::setup() {
    /*
     * \brief uniform values are treated immediatly. For spatially variable
     * fields, we return the information needed to evaluate them
     */
    // This lambda function builds an *evaluator*. Its role is to fill up
    // the `v` vector
    auto dispatch =
        [](std::vector<real>& v,
           std::map<std::string,
                    mgis::variant<real, mgis::span<real>, std::vector<real>>>&
               values,
           const std::vector<mgis::behaviour::Variable>& ds) {
          mgis::raise_if(ds.size() != v.size(),
                         "integrate: ill allocated memory");
          // evaluators
          std::vector<std::tuple<size_type, real*>> evs;
          auto i = mgis::size_type{};
          for (const auto& d : ds) {
            if (d.type != mgis::behaviour::Variable::SCALAR) {
              mgis::raise("integrate: invalid type for variable '" + d.name +
                          "'");
            }
            auto p = values.find(d.name);
            if (p == values.end()) {
              auto msg = std::string{"integrate: no variable named '" + d.name +
                                     "' declared"};
              if (!values.empty()) {
                msg += "\nThe following variables were declared: ";
                for (const auto& vs : values) {
                  msg += "\n- " + vs.first;
                }
              } else {
                msg += "\nNo variable declared.";
              }
              mgis::raise(msg);
            }
	    // depending on the type of p->second, we are branching
	    // on one of the following procedure:
            if (mgis::holds_alternative<real>(p->second)) {
	      // if uniform field, copy p->second into v[i]
	      // `evs` will be empty.
              v[i] = mgis::get<real>(p->second);
            } else if (mgis::holds_alternative<mgis::span<real>>(p->second)) {
	      // if we have a span, we store in evs this span for future use
              evs.push_back(std::make_tuple(
                  i, mgis::get<mgis::span<real>>(p->second).data()));
            } else {
	      // if we have a vector, we store in evs this vector for future use
              evs.push_back(std::make_tuple(
                  i, mgis::get<std::vector<real>>(p->second).data()));
            }
            ++i;
          }
          return evs;
        };  // end of dispatch

    // The `b` field refers to MGIS MaterialDataManager class, which
    // is an ancestor class of the current BehaviourIntegratorBase class.
    this->wks.mps_evaluators =
        dispatch(this->wks.mps, this->s1.material_properties, this->b.mps);
    this->wks.esvs0_evaluators = dispatch(
        this->wks.esvs0, this->s0.external_state_variables, this->b.esvs);
    this->wks.esvs1_evaluators = dispatch(
        this->wks.esvs1, this->s1.external_state_variables, this->b.esvs);
  }  // end of revert

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
    const auto g_offset = this->s0.gradients_stride * ip;
    const auto t_offset = this->s0.thermodynamic_forces_stride * ip;
    const auto isvs_offset = this->s0.internal_state_variables_stride * ip;
    // Fill vector `v` with material properties or external state variables.
    // Evaluators are used to deal with both uniform or spatially
    // variable quantities.
    auto eval = [](std::vector<real>& v,
                   const std::vector<std::tuple<size_type, real*>>& evs,
                   const size_type i) {
      for (const auto& ev : evs) {
        v[std::get<0>(ev)] = std::get<1>(ev)[i];
      }
    };  // end of eval
    eval(this->wks.mps, this->wks.mps_evaluators, ip);
    eval(this->wks.esvs0, this->wks.esvs0_evaluators, ip);
    eval(this->wks.esvs1, this->wks.esvs1_evaluators, ip);
    //
    mgis::behaviour::BehaviourDataView v;
    v.rdt = real(1);
    v.dt = this->time_increment;
    v.K = this->K.data() + this->K_stride * ip;
    v.s0.gradients = this->s0.gradients.data() + g_offset;
    v.s1.gradients = this->s1.gradients.data() + g_offset;
    v.s0.thermodynamic_forces =
        this->s0.thermodynamic_forces.data() + t_offset;
    v.s1.thermodynamic_forces =
        this->s1.thermodynamic_forces.data() + t_offset;
    v.s0.material_properties = this->wks.mps.data();
    v.s1.material_properties = this->wks.mps.data();
    v.s0.internal_state_variables =
        this->s0.internal_state_variables.data() + isvs_offset;
    v.s1.internal_state_variables =
        this->s1.internal_state_variables.data() + isvs_offset;
    if (this->b.computesStoredEnergy) {
      v.s0.stored_energy = this->s0.stored_energies.data() + ip;
      v.s1.stored_energy = this->s1.stored_energies.data() + ip;
    } else {
      v.s0.stored_energy = nullptr;
      v.s1.stored_energy = nullptr;
    }
    if (this->b.computesDissipatedEnergy) {
      v.s0.dissipated_energy = this->s0.dissipated_energies.data() + ip;
      v.s1.dissipated_energy = this->s1.dissipated_energies.data() + ip;
    } else {
      v.s0.dissipated_energy = nullptr;
      v.s1.dissipated_energy = nullptr;
    }
    v.s0.external_state_variables = this->wks.esvs0.data();
    v.s1.external_state_variables = this->wks.esvs1.data();
    v.K[0] = static_cast<int>(it);
    const auto r = mgis::behaviour::integrate(v, this->b);
    if ((r != 0) && (r != 1)) {
      mgis::raise(
          "BehaviourIntegratorBase::integrate: "
          "behaviour integration failed (" +
          std::to_string(r) + ")");
    }
  }  // end of BehaviourIntegratorBase::integrate

  void BehaviourIntegratorBase::revert(){
    mgis::behaviour::revert(*this);
  }  // end of revert

  void BehaviourIntegratorBase::update(){
    mgis::behaviour::update(*this);
  }  // end of update

  BehaviourIntegratorBase::~BehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis
