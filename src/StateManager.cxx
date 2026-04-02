/*!
 * \file   StateManager.cxx
 * \brief  This file implements the `StateManager` class
 * \author Thomas Helfer
 * \date   31/03/2026
 */

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/StateManager.hxx"

namespace mfem_mgis {

  bool StateManager::PartialQuadratureFunctionManager::add(
      Context& ctx,
      std::string_view n,
      ImmutablePartialQuadratureFunctionView f,
      const TimeStepStage ts) noexcept {
    auto& m = ts == bts ? this->qfunctions_bts : this->qfunctions_ets;
    if (m.find(n) != m.end()) {
      return ctx.registerErrorMessage("a partial quadrature function named '" +
                                      std::string{n} +
                                      "' has already been registred");
    }
    m.insert({std::string{n}, f});
    return true;
  }  // end of add

  std::optional<ImmutablePartialQuadratureFunctionView>
  StateManager::PartialQuadratureFunctionManager::get(
      Context& ctx, std::string_view n, const TimeStepStage ts) const noexcept {
    const auto& m = ts == bts ? this->qfunctions_bts : this->qfunctions_ets;
    const auto p = m.find(n);
    if (p == m.end()) {
      auto msg = std::string{"no partial quadrature function named '" +
                             std::string{n} + "' registred"};
      if (!m.empty()) {
        msg += ". The following partial quadrature functions are declared:";
        for (const auto& [name, f] : m) {
          static_cast<void>(f);
          msg += "\n- '" + name + "'";
        }
      }
      return ctx.registerErrorMessage(msg);
    }
    return p->second;
  }  // end of get

  bool StateManager::PartialQuadratureFunctionManager::contains(
      std::string_view n, const TimeStepStage ts) const noexcept {
    const auto& m = ts == bts ? this->qfunctions_bts : this->qfunctions_ets;
    return m.find(n) != m.end();
  }  // end of get

  StateManager::StateManager(const MeshDiscretization& m) noexcept
      : PartialQuadratureSpaceIdentifiersManager(m) {}

  bool StateManager::add(Context& ctx,
                         std::string_view n,
                         const std::shared_ptr<PartialQuadratureFunction>& f,
                         const TimeStepStage ts) noexcept {
    if (f.get() == nullptr) {
      return ctx.registerErrorMessage("invalid partial quadrature function '" +
                                      std::string{n} + "'");
    }
    const auto qspace = f->getPartialQuadratureSpacePointer();
    const auto success = this->add(ctx, n, *f, ts);
    if (success) {
      this->qfunctions_pointers.push_back(f);
    }
    return success;
  }  // end of add

  bool StateManager::add(Context& ctx,
                         std::string_view n,
                         ImmutablePartialQuadratureFunctionView f,
                         const TimeStepStage ts) noexcept {
    const auto qspace = f.getPartialQuadratureSpacePointer();
    const auto oqid = this->getIdentifier(ctx, qspace);
    if (isInvalid(oqid)) {
      return false;
    }
    const auto id = std::make_pair(qspace->getId(), *oqid);
    auto p = this->qfunctions.find(id);
    if (p == this->qfunctions.end()) {
      auto ptr = std::make_unique<PartialQuadratureFunctionManager>();
      p = this->qfunctions.insert({id, std::move(ptr)}).first;
    }
    return p->second->add(ctx, n, f, ts);
  }  // end of add

  std::optional<bool> StateManager::contains(
      Context& ctx,
      const std::shared_ptr<const PartialQuadratureSpace>& qspace,
      std::string_view n,
      const TimeStepStage ts) const noexcept {
    const auto oqid =
        static_cast<const StateManager&>(*this).getIdentifier(ctx, qspace);
    if (isInvalid(oqid)) {
      return {};
    }
    const auto id = std::make_pair(qspace->getId(), *oqid);
    auto p = this->qfunctions.find(id);
    if (p == this->qfunctions.end()) {
      return false;
    }
    return p->second->contains(n, ts);
  }  // end of get

  std::optional<ImmutablePartialQuadratureFunctionView> StateManager::get(
      Context& ctx,
      const std::shared_ptr<const PartialQuadratureSpace>& qspace,
      std::string_view n,
      const TimeStepStage ts) const noexcept {
    const auto oqid =
        static_cast<const StateManager&>(*this).getIdentifier(ctx, qspace);
    if (isInvalid(oqid)) {
      return {};
    }
    const auto id = std::make_pair(qspace->getId(), *oqid);
    auto p = this->qfunctions.find(id);
    if (p == this->qfunctions.end()) {
      return ctx.registerErrorMessage("no partial function named '" +
                                      std::string{n} + "' on material '" +
                                      std::to_string(qspace->getId()) + "'");
    }
    return p->second->get(ctx, n, ts);
  }  // end of get

  StateManager::~StateManager() noexcept = default;

  static bool addPartialQuadratureFunctions(
      Context& ctx, StateManager& s, const AbstractBehaviourIntegrator& bi) {
    if (!bi.hasMaterial()) {
      return true;
    }
    const auto om = bi.getMaterial(ctx);
    if (isInvalid(om)) {
      return false;
    }
    const auto qspace = om->getPartialQuadratureSpacePointer();
    const auto& b = om->b;  // underlying behaviour
    // declaring gradients
    for (const auto& g : b.gradients) {
      const auto og_bts = getGradient(ctx, *om, g.name, bts);
      const auto og_ets = getGradient(ctx, *om, g.name, ets);
      if (!areValid(og_bts, og_ets)) {
        return false;
      }
      if (!s.add(ctx, g.name, *og_bts, bts)) {
        return false;
      }
      if (!s.add(ctx, g.name, *og_ets, ets)) {
        return false;
      }
    }
    // declaring thermodynamics forces
    for (const auto& th : b.thermodynamic_forces) {
      const auto oth_bts = getThermodynamicForce(ctx, *om, th.name, bts);
      const auto oth_ets = getThermodynamicForce(ctx, *om, th.name, ets);
      if (!areValid(oth_bts, oth_ets)) {
        return false;
      }
      if (!s.add(ctx, th.name, *oth_bts, bts)) {
        return false;
      }
      if (!s.add(ctx, th.name, *oth_ets, ets)) {
        return false;
      }
    }
    // declaring internal state variables
    for (const auto& isv : b.isvs) {
      const auto oisv_bts = getInternalStateVariable(ctx, *om, isv.name, bts);
      const auto oisv_ets = getInternalStateVariable(ctx, *om, isv.name, ets);
      if (!areValid(oisv_bts, oisv_ets)) {
        return false;
      }
      if (!s.add(ctx, isv.name, *oisv_bts, bts)) {
        return false;
      }
      if (!s.add(ctx, isv.name, *oisv_ets, ets)) {
        return false;
      }
    }
    //
    return true;
  }  // end of addPartialQuadratureFunctions

  bool addPartialQuadratureFunctions(
      Context& ctx,
      StateManager& s,
      const AbstractNonLinearEvolutionProblem& p) noexcept {
    // loop on materials for which at least on behaviour integrator is defined
    for (const auto m : p.getAssignedMaterialsIdentifiers()) {
      const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, m);
      if (isInvalid(onbis)) {
        return false;
      }
      for (size_type n = 0; n != onbis; ++n) {
        const auto obi = p.getBehaviourIntegrator(ctx, m, n);
        if (isInvalid(obi)) {
          return false;
        }
        if (!addPartialQuadratureFunctions(ctx, s, *obi)) {
          return false;
        }
      }
    }
    return true;
  }  // end of addPartialQuadratureFunctions

}  // end of namespace mfem_mgis