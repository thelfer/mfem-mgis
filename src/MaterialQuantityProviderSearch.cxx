/*!
 * \file   src/MaterialQuantityProviderSearch.cxx
 * \brief
 * \date   23/09/2026
 */

#include <vector>
#include "MGIS/Behaviour/Variable.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/MaterialQuantityProviderSearch.hxx"

namespace mfem_mgis {

  /*!
   * \brief search the behaviour integrator providing a quantity on a location
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] l: location
   * \param[in] n: name of the quantity
   * \param[in] get_variables: function returning the list of variables of a
   * material in which the quantity is searched
   */
  template <typename VariablesGetter>
  [[nodiscard]] static std::optional<MaterialQuantityProviderSearchResult>
  searchMaterialQuantityProvider(
      Context& ctx,
      const AbstractNonLinearEvolutionProblem& p,
      const LocationIdentifier& l,
      std::string_view n,
      const VariablesGetter& get_variables) noexcept {
    using Status = MaterialQuantityProviderSearchResult::Status;
    if (isInvalid(l)) {
      return ctx.registerErrorMessage("invalid location identifier");
    }
    auto r = MaterialQuantityProviderSearchResult{
        .behaviour_integrator = {}, .status = Status::NO_PROVIDER};
    if (!l.material_identifier.has_value()) {
      // behaviour integrators are currently only defined on materials
      return r;
    }
    const auto mid = l.material_identifier->id;
    const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, mid);
    if (isInvalid(onbis)) {
      return {};
    }
    for (auto i = size_type{}; i != *onbis; ++i) {
      const auto obi = p.getBehaviourIntegrator(ctx, mid, i);
      if (isInvalid(obi)) {
        return {};
      }
      if (!obi->hasMaterial()) {
        continue;
      }
      const auto om = obi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      if (!contains(get_variables(*om), n)) {
        continue;
      }
      if (r.status == Status::SUCCESS) {
        return MaterialQuantityProviderSearchResult{
            .behaviour_integrator = {}, .status = Status::MULTIPLE_PROVIDERS};
      }
      r.behaviour_integrator = obi;
      r.status = Status::SUCCESS;
    }
    return r;
  }  // end of searchMaterialQuantityProvider

  std::optional<MaterialQuantityProviderSearchResult> hasGradientProvider(
      Context& ctx,
      const AbstractNonLinearEvolutionProblem& p,
      const LocationIdentifier& l,
      std::string_view n) noexcept {
    return searchMaterialQuantityProvider(
        ctx, p, l, n,
        [](const Material& m) -> const std::vector<mgis::behaviour::Variable>& {
          return m.b.gradients;
        });
  }  // end of hasGradientProvider

  std::optional<MaterialQuantityProviderSearchResult>
  hasThermodynamicForceProvider(Context& ctx,
                                const AbstractNonLinearEvolutionProblem& p,
                                const LocationIdentifier& l,
                                std::string_view n) noexcept {
    return searchMaterialQuantityProvider(
        ctx, p, l, n,
        [](const Material& m) -> const std::vector<mgis::behaviour::Variable>& {
          return m.b.thermodynamic_forces;
        });
  }  // end of hasThermodynamicForceProvider

  std::optional<MaterialQuantityProviderSearchResult>
  hasInternalStateVariableProvider(Context& ctx,
                                   const AbstractNonLinearEvolutionProblem& p,
                                   const LocationIdentifier& l,
                                   std::string_view n) noexcept {
    return searchMaterialQuantityProvider(
        ctx, p, l, n,
        [](const Material& m) -> const std::vector<mgis::behaviour::Variable>& {
          return m.b.isvs;
        });
  }  // end of hasInternalStateVariableProvider

}  // end of namespace mfem_mgis
