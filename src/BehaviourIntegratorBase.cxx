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
#include "MFEMMGIS/AbstractPartialQuadratureFunctionEvaluator.hxx"
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

  bool BehaviourIntegratorBase::setMaterialProperty(
      Context& ctx,
      std::string_view name,
      std::shared_ptr<const AbstractPartialQuadratureFunctionEvaluator> e,
      const TimeStepStage ts) noexcept {
    const auto omp = getVariable(ctx, this->b.mps, name);
    if (isInvalid(omp)) {
      return false;
    }
    if (e.get() == nullptr) {
      return ctx.registerErrorMessage(
          "invalid partial quadrature function evaluator given for material "
          "property '" +
          (*omp)->name + "'");
    }
    if ((*omp)->type != mgis::behaviour::Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "invalid material property '" + (*omp)->name +
          "' (only scalar material properties are supported)");
    }
    if (e->getNumberOfComponents() != 1) {
      return ctx.registerErrorMessage(
          "invalid number of components for partial quadrature function "
          "evaluator given for material property '" +
          (*omp)->name + "' (" + std::to_string(e->getNumberOfComponents()) +
          ", expected 1)");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      if (!unsetMaterialProperty(ctx, this->s0, name)) {
        return false;
      }
      this->material_properties_evaluators_bts.insert_or_assign(
          std::string{name}, e);
    } else {
      if (!unsetMaterialProperty(ctx, this->s1, name)) {
        return false;
      }
      this->material_properties_evaluators_ets.insert_or_assign(
          std::string{name}, e);
    }
    return true;
  }  // end of setMaterialProperty

  bool BehaviourIntegratorBase::setExternalStateVariable(
      Context& ctx,
      std::string_view name,
      std::shared_ptr<const AbstractPartialQuadratureFunctionEvaluator> e,
      const TimeStepStage ts) noexcept {
    const auto oesv = getVariable(ctx, this->b.esvs, name);
    if (isInvalid(oesv)) {
      return false;
    }
    if (e.get() == nullptr) {
      return ctx.registerErrorMessage(
          "invalid partial quadrature function evaluator given for external "
          "state variable '" +
          (*oesv)->name + "'");
    }
    const auto osize = getVariableSize(ctx, *(*oesv), this->b.hypothesis);
    if (isInvalid(osize)) {
      return false;
    }
    if (e->getNumberOfComponents() != *osize) {
      return ctx.registerErrorMessage(
          "invalid number of components for partial quadrature function "
          "evaluator given for external state variable '" +
          (*oesv)->name + "' (" + std::to_string(e->getNumberOfComponents()) +
          ", expected " + std::to_string(*osize) + ")");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      if (!unsetExternalStateVariable(ctx, this->s0, name)) {
        return false;
      }
      this->external_state_variables_evaluators_bts.insert_or_assign(
          std::string{name}, e);
    } else {
      if (!unsetExternalStateVariable(ctx, this->s1, name)) {
        return false;
      }
      this->external_state_variables_evaluators_ets.insert_or_assign(
          std::string{name}, e);
    }
    return true;
  }  // end of setExternalStateVariable

  void BehaviourIntegratorBase::throwInvalidBehaviourType(
      const char* const mn, const char* const m) const {
    auto msg = std::string(mn) + ": invalid behaviour type";
    if (m != nullptr) {
      msg += "(" + std::string(m) + ')';
    }
    raise(msg);
  }  // end of throwInvalidBehaviourType

  void BehaviourIntegratorBase::throwInvalidBehaviourKinematic(
      const char* const mn, const char* const m) const {
    auto msg = std::string(mn) + ": invalid behaviour type";
    if (m != nullptr) {
      msg += "(" + std::string(m) + ')';
    }
    raise(msg);
  }  // end of throwInvalidBehaviourType

  const PartialQuadratureSpace&
  BehaviourIntegratorBase::getPartialQuadratureSpace() const noexcept {
    return this->getMaterial().getPartialQuadratureSpace();
  }  // end of getPartialQuadratureSpace

  real BehaviourIntegratorBase::getTimeIncrement() const noexcept {
    return this->time_increment;
  }  // end of getTimeIncrement

  void BehaviourIntegratorBase::setTimeIncrement(const real dt) {
    this->time_increment = dt;
  }  // end of setTimeIncrement

  bool BehaviourIntegratorBase::hasMaterial() const noexcept { return true; }

  OptionalReference<Material> BehaviourIntegratorBase::getMaterial(
      Context&) noexcept {
    return OptionalReference<Material>{this};
  }  // end of getMaterial

  OptionalReference<const Material> BehaviourIntegratorBase::getMaterial(
      Context&) const noexcept {
    return OptionalReference<const Material>{this};
  }  // end of getMaterial

  Material& BehaviourIntegratorBase::getMaterial() {
    return *this;
  }  // end of getMaterial

  const Material& BehaviourIntegratorBase::getMaterial() const {
    return *this;
  }  // end of getMaterial

  bool BehaviourIntegratorBase::setup(Context& ctx,
                                      const real,
                                      const real) noexcept {
    this->wks.pqfcts_mps_bts.clear();
    this->wks.pqfcts_mps_ets.clear();
    this->wks.pqfcts_esvs_bts.clear();
    this->wks.pqfcts_esvs_ets.clear();
    for (const auto& [name, e] : this->material_properties_evaluators_bts) {
      auto f = evaluate(ctx, *e);
      if (isInvalid(f)) {
        return false;
      }
      this->wks.pqfcts_mps_bts.insert_or_assign(name, f);
      if (!mgis::behaviour::setMaterialProperty(
              ctx, this->s0, name, f->getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->material_properties_evaluators_ets) {
      auto f = evaluate(ctx, *e);
      if (isInvalid(f)) {
        return false;
      }
      this->wks.pqfcts_mps_ets.insert_or_assign(name, f);
      if (!mgis::behaviour::setMaterialProperty(
              ctx, this->s1, name, f->getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->external_state_variables_evaluators_bts) {
      auto f = evaluate(ctx, *e);
      if (isInvalid(f)) {
        return false;
      }
      this->wks.pqfcts_mps_bts.insert_or_assign(name, f);
      if (!mgis::behaviour::setExternalStateVariable(
              ctx, this->s0, name, f->getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->external_state_variables_evaluators_ets) {
      auto f = evaluate(ctx, *e);
      if (isInvalid(f)) {
        return false;
      }
      this->wks.pqfcts_mps_ets.insert_or_assign(name, f);
      if (!mgis::behaviour::setExternalStateVariable(
              ctx, this->s1, name, f->getValues(),
              mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
              mgis::behaviour::MaterialStateManager::NOUPDATE)) {
        return false;
      }
    }
    /*
     * \brief uniform values are treated immediatly. For spatially variable
     * fields, we return the information needed to evaluate them
     */
    // This lambda function builds an *evaluator*. Its role is to fill up
    // the `v` vector
    auto dispatch =
        [this, &ctx](
            std::vector<real>& v,
            const std::map<std::string,
                           mgis::behaviour::MaterialStateManager::FieldHolder,
                           std::less<>>& field_holders,
            const std::vector<mgis::behaviour::Variable>& ds)
        -> std::optional<
            std::vector<std::tuple<size_type, size_type, const real*>>> {
      const auto h = this->getMaterial().b.hypothesis;
      raise_if(v.size() != mgis::behaviour::getArraySize(ds, h),
               "integrate: ill allocated memory");
      // evaluators
      auto evs = std::vector<std::tuple<size_type, size_type, const real*>>{};
      auto i = mgis::size_type{};
      for (const auto& d : ds) {
        const auto vsize = mgis::behaviour::getVariableSize(d, h);
        auto p = field_holders.find(d.name);
        if (p == field_holders.end()) {
          auto msg = std::string{"integrate: no variable named '" + d.name +
                                 "' declared"};
          if (!field_holders.empty()) {
            msg += "\nThe following variables were declared: ";
            for (const auto& vs : field_holders) {
              msg += "\n- " + vs.first;
            }
          } else {
            msg += "\nNo variable declared.";
          }
          return ctx.registerErrorMessage(msg);
        }
        // depending on the type of p->second, we are branching
        // on one of the following procedure:
        const auto& field_value = p->second.value;
        if (std::holds_alternative<real>(field_value)) {
          if (vsize != 1) {
            return ctx.registerErrorMessage(
                "invalid number of values given for variable '" + d.name + "'");
          }
          // if uniform field, copy field_value into v[i]
          // `evs` will be untouched.
          v[i] = std::get<real>(field_value);
        } else if (std::holds_alternative<std::span<real>>(field_value)) {
          const auto& variable_values = std::get<std::span<real>>(field_value);
          if (variable_values.size() == v.size()) {
            std::copy(variable_values.begin(), variable_values.end(),
                      v.begin() + i);
          } else {
            if (variable_values.size() != vsize * (this->n)) {
              return ctx.registerErrorMessage(
                     "invalid number of variable values given for variable '" +
                         d.name + "'");
            }
            // if we have a span, we store in evs this span for future use
            evs.push_back(std::make_tuple(i, vsize, variable_values.data()));
          }
        } else {
          const auto& variable_values =
              std::get<std::vector<real>>(field_value);
          if (variable_values.size() == v.size()) {
            std::copy(variable_values.begin(), variable_values.end(),
                      v.begin() + i);
          } else {
            if(variable_values.size() != vsize * (this->n)){
              return ctx.registerErrorMessage(
                  "invalid number of variable values given for variable '" +
                  d.name + "'");
            }
            // if we have a vector, we store in evs this vector for future
            // use
            evs.push_back(std::make_tuple(i, vsize, variable_values.data()));
          }
        }
        i += vsize;
      }
      return evs;
    };  // end of dispatch
    // The `b` field refers to MGIS MaterialDataManager class, which
    // is an ancestor class of the current BehaviourIntegratorBase class.
    auto oevaluators =
        dispatch(this->wks.mps, this->s1.material_properties, this->b.mps);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.mps_evaluators = *oevaluators;
    // external state variables at the beginning of the time step
    oevaluators = dispatch(this->wks.esvs0, this->s0.external_state_variables,
                           this->b.esvs);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.esvs0_evaluators = *oevaluators;
    // external state variables at the end of the time step
    oevaluators = dispatch(this->wks.esvs1, this->s1.external_state_variables,
                           this->b.esvs);
    if (isInvalid(oevaluators)) {
      return false;
    }
    this->wks.esvs1_evaluators = *oevaluators;
    return true;
  }  // end of setup

  void BehaviourIntegratorBase::checkHypotheses(const Hypothesis h) const {
    using namespace mgis::behaviour;
    if (this->b.hypothesis != h) {
      const auto h1 = std::string(toString(this->b.hypothesis));
      const auto h2 = std::string(toString(h));
      raise(
          "BehaviourIntegratorBase::checkHypotheses: "
          "the behaviour hypothesis (" +
          h1 + ") does not match the integrator hypothesis (" + h2 + ")");
    }
  }  // end of BehaviourIntegratorBase::checkHypotheses

  bool BehaviourIntegratorBase::performsLocalBehaviourIntegration(
      const size_type ip, const IntegrationType it) {
    char error_msg[512];
    error_msg[0] = '\0';
    const auto g_offset = this->s0.gradients_stride * ip;
    const auto t_offset = this->s0.thermodynamic_forces_stride * ip;
    const auto isvs_offset = this->s0.internal_state_variables_stride * ip;
    // Fill vector `v` with material properties or external state variables.
    // Evaluators are used to deal with both uniform or spatially
    // variable quantities.
    auto eval =
        [](std::vector<real>& v,
           const std::vector<std::tuple<size_type, size_type, const real*>>&
               evs,
           const size_type i) {
          for (const auto& [offset, size, values] : evs) {
            if (size == 1) {
              v[offset] = values[i];
            } else {
              const auto b_ptr = values + i * size;
              const auto e_ptr = values + (i + 1) * size;
              std::copy(b_ptr, e_ptr, v.begin() + offset);
            }
          }
        };  // end of eval
    eval(this->wks.mps, this->wks.mps_evaluators, ip);
    eval(this->wks.esvs0, this->wks.esvs0_evaluators, ip);
    eval(this->wks.esvs1, this->wks.esvs1_evaluators, ip);
    //
    rdt = real{1};
    mgis::behaviour::BehaviourDataView v;
    v.error_message = error_msg;
    v.rdt = &rdt;
    v.dt = this->time_increment;
    v.K = this->K.data() + this->K_stride * ip;
    v.speed_of_sound = nullptr;
    v.s0.gradients = this->s0.gradients.data() + g_offset;
    v.s1.gradients = this->s1.gradients.data() + g_offset;
    v.s0.thermodynamic_forces = this->s0.thermodynamic_forces.data() + t_offset;
    v.s1.thermodynamic_forces = this->s1.thermodynamic_forces.data() + t_offset;
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
    v.s0.mass_density = nullptr;
    v.s1.mass_density = nullptr;
    v.s0.external_state_variables = this->wks.esvs0.data();
    v.s1.external_state_variables = this->wks.esvs1.data();
    v.K[0] = static_cast<int>(it);
    const auto r = mgis::behaviour::integrate(v, this->b);
    if (!((r == 0) || (r == 1))) {
      std::cerr << "behaviour integration failed: " << error_msg
                << "(exit value: " << r << ")" << std::endl;
    }
    return (r == 0) || (r == 1);
  }  // end of BehaviourIntegratorBase::integrate

  bool BehaviourIntegratorBase::cleanup(Context& ctx) noexcept {
    for (const auto& [name, e] : this->material_properties_evaluators_bts) {
      std::ignore = e;
      if (!mgis::behaviour::unsetMaterialProperty(ctx, this->s0, name)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->material_properties_evaluators_ets) {
      std::ignore = e;
      if (!mgis::behaviour::unsetMaterialProperty(
              ctx, this->s1, name)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->external_state_variables_evaluators_bts) {
      std::ignore = e;
      if (!mgis::behaviour::unsetExternalStateVariable(ctx, this->s0, name)) {
        return false;
      }
    }
    for (const auto& [name, e] : this->external_state_variables_evaluators_ets) {
      std::ignore = e;
      if (!mgis::behaviour::unsetExternalStateVariable(ctx, this->s1, name)) {
        return false;
      }
    }
    this->wks.pqfcts_mps_bts.clear();
    this->wks.pqfcts_mps_ets.clear();
    this->wks.pqfcts_esvs_bts.clear();
    this->wks.pqfcts_esvs_ets.clear();
    return true;
  } // end of cleanup

  void BehaviourIntegratorBase::revert() {
    mgis::behaviour::revert(*this);
  }  // end of revert

  void BehaviourIntegratorBase::update() {
    mgis::behaviour::update(*this);
  }  // end of update

  void BehaviourIntegratorBase::setMacroscopicGradients(
      std::span<const real> g) {
    Material::setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  bool BehaviourIntegratorBase::requiresCurrentSolutionForResidualAssembly()
      const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForResidualAssembly

  bool BehaviourIntegratorBase::requiresCurrentSolutionForJacobianAssembly()
      const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForJacobianAssembly

  BehaviourIntegratorBase::~BehaviourIntegratorBase() = default;

}  // end of namespace mfem_mgis
