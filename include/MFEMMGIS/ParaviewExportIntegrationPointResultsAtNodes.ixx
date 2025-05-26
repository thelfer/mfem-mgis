/*!
 * \file   include/MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.ixx
 * \brief
 * \author Thomas Helfer
 * \date   18/08/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        cycle(0) {
    std::cerr << "ParaviewExportIntegrationPointResultsAtNodes\n";
    checkParameters(params, {"OutputFileName", "Materials", "Results"});
    // if Materials exists, use it, otherwise, take all materials
    this->materials_identifiers = getMaterialsIdentifiers(p, params);
    if (!contains(params, "Results")) {
      raise(
          "ParaviewExportIntegrationPointResultsAtNodes::"
          "ParaviewExportIntegrationPointResultsAtNodes: "
          "no results to export declared");
    }
    //
    auto add_result = [this, &p](const std::string& rn) {
      std::cerr << "tata" << rn << "\n";
      auto r = IntegrationPointResult{};
      r.name = rn;
      this->getResultDescription(r, p);
      std::tie(r.fespace, r.f) = makeGridFunction<parallel>(
          this->getPartialQuadratureFunctionViews(p, r));
      // registring
      this->exporter.RegisterField(r.name, r.f.get());
      // saving
      this->results.push_back(std::move(r));
    };
    if (is<std::string>(params, "Results")) {
      add_result(get<std::string>(params, "Results"));
    } else {
      for (const auto& rn : get<std::vector<Parameter>>(params, "Results")) {
        add_result(get<std::string>(rn));
      }
    }
  }  // end of ParaviewExportIntegrationPointResultsAtNodes

  template <bool parallel>
  void ParaviewExportIntegrationPointResultsAtNodes<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);
    // updating grid functions
    for (auto& r : this->results) {
      updateGridFunction<parallel>(
          *(r.f), this->getPartialQuadratureFunctionViews(p, r));
    }
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  void
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::getResultDescription(
      IntegrationPointResult& r,
      const NonLinearEvolutionProblemImplementation<parallel>& p) {
    using namespace mgis::behaviour;
    auto c = typename IntegrationPointResult::Category{};
    auto s = size_type{};
    auto first = true;
    if (this->materials_identifiers.empty()) {
      raise("getResultDescription: empty list of material identifiers");
    }
    for (const auto& mid : this->materials_identifiers) {
      const auto [c2, s2] = [
        &p, &mid, &r
      ]() -> std::tuple<typename IntegrationPointResult::Category, size_type> {
        const auto& m = p.getMaterial(mid);
        const auto& h = m.b.hypothesis;
        if (contains(m.b.gradients, r.name)) {
          return {IntegrationPointResult::GRADIENTS,
                  getVariableSize(getVariable(m.b.gradients, r.name), h)};
        } else if (contains(m.b.thermodynamic_forces, r.name)) {
          return {IntegrationPointResult::THERMODYNAMIC_FORCES,
                  getVariableSize(getVariable(m.b.thermodynamic_forces, r.name),
                                  h)};
        } else if (!contains(m.b.isvs, r.name)) {
          raise("getResultDescription: no result '" + std::string(r.name) +
                "' found for material '" + std::to_string(mid) + "'");
        }
        return {IntegrationPointResult::INTERNAL_STATE_VARIABLES,
                getVariableSize(getVariable(m.b.isvs, r.name), h)};
      }
      ();
      if (first) {
        s = s2;
        c = c2;
      } else {
        if (s != s2) {
          raise(
              "getResultDescription: inconsistent number of components of "
              "field '" +
              std::string(r.name) + "'");
        }
        if (c != c2) {
          raise(
              "getResultDescription: inconsistent nature for "
              "field '" +
              std::string(r.name) + "'");
        }
      }
    }
    r.category = c;
    r.number_of_components = s;
  }  // end of getResultDescription

  template <bool parallel>
  std::vector<ImmutablePartialQuadratureFunctionView>
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      getPartialQuadratureFunctionViews(
          const NonLinearEvolutionProblemImplementation<parallel>& p,
          const IntegrationPointResult& r) {
    auto fcts = std::vector<ImmutablePartialQuadratureFunctionView>{};
    // creating the partial functions per materials
    for (const auto& mid : this->materials_identifiers) {
      const auto& m = p.getMaterial(mid);
      if (r.category == IntegrationPointResult::GRADIENTS) {
        fcts.push_back(getGradient(m, r.name));
      } else if (r.category == IntegrationPointResult::THERMODYNAMIC_FORCES) {
        fcts.push_back(getThermodynamicForce(m, r.name));
      } else {
        fcts.push_back(getInternalStateVariable(m, r.name));
      }
    }
    return fcts;
  }  // end of getPartialQuadratureFunctionViews

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodes<
      parallel>::~ParaviewExportIntegrationPointResultsAtNodes() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX */
