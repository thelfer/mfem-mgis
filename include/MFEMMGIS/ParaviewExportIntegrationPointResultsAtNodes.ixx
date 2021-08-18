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
  struct ParaviewExportIntegrationPointResultsAtNodes<
      parallel>::ResultCoefficientBase {
    //
    ResultCoefficientBase(const IntegrationPointResult& r,
                          NonLinearEvolutionProblemImplementation<parallel>& p,
                          const std::vector<size_type>& mids) {
      // creating the partial functions per materials
      for (const auto& mid : mids) {
        auto& m = p.getMaterial(mid);
        if (r.category == IntegrationPointResult::GRADIENTS) {
          this->functions.insert({mid, getGradient(m, r.name)});
        } else if (r.category == IntegrationPointResult::THERMODYNAMIC_FORCES) {
          this->functions.insert({mid, getThermodynamicForce(m, r.name)});
        } else {
          this->functions.insert({mid, getInternalStateVariable(m, r.name)});
        }
      }
    }
    //
    ResultCoefficientBase(ResultCoefficientBase&&) = default;
    ResultCoefficientBase(const ResultCoefficientBase&) = default;
    ResultCoefficientBase& operator=(ResultCoefficientBase&&) = default;
    ResultCoefficientBase& operator=(const ResultCoefficientBase&) = default;

   protected:
    std::map<size_type, PartialQuadratureFunction> functions;
  };  // end of ResultCoefficientBase

  template <bool parallel>
  struct ParaviewExportIntegrationPointResultsAtNodes<
      parallel>::ScalarResultCoefficient : public ResultCoefficientBase,
                                           public mfem::Coefficient {
    //
    ScalarResultCoefficient(
        const IntegrationPointResult& r,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const std::vector<size_type>& mids)
        : ResultCoefficientBase(r, p, mids) {}
    //
    ScalarResultCoefficient(ScalarResultCoefficient&&) = default;
    ScalarResultCoefficient(const ScalarResultCoefficient&) = default;
    ScalarResultCoefficient& operator=(ScalarResultCoefficient&&) = default;
    ScalarResultCoefficient& operator=(const ScalarResultCoefficient&) =
        default;
    //
    double Eval(mfem::ElementTransformation& tr,
                const mfem::IntegrationPoint& i) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        return 0;
      }
      return p->second.getIntegrationPointValue(tr.ElementNo, i.index);
    }  // end of Eval
  };

  template <bool parallel>
  struct ParaviewExportIntegrationPointResultsAtNodes<
      parallel>::MultiComponentsResultCoefficient
      : public ResultCoefficientBase,
        public mfem::VectorCoefficient {
    //
    MultiComponentsResultCoefficient(
        const IntegrationPointResult& r,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const std::vector<size_type>& mids)
        : ResultCoefficientBase(r, p, mids),
          mfem::VectorCoefficient(r.number_of_components) {}
    //
    MultiComponentsResultCoefficient(MultiComponentsResultCoefficient&&) =
        default;
    MultiComponentsResultCoefficient(const MultiComponentsResultCoefficient&) =
        default;
    MultiComponentsResultCoefficient& operator=(
        MultiComponentsResultCoefficient&&) = default;
    MultiComponentsResultCoefficient& operator=(
        const MultiComponentsResultCoefficient&) = default;
    //
    void Eval(mfem::Vector& values,
              mfem::ElementTransformation& tr,
              const mfem::IntegrationPoint& ip) override {
      const auto mid = tr.Attribute;
      const auto p = this->functions.find(mid);
      if (p == this->functions.end()) {
        values = 0.;
      } else {
        const auto rvalues =
            p->second.getIntegrationPointValues(tr.ElementNo, ip.index);
        for (size_type i = 0; i != this->GetVDim(); ++i) {
          values[i] = rvalues[i];
        }
      }
    }  // end of Eval
  };

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        cycle(0) {
    checkParameters(params, {"OutputFileName", "Materials", "Results"});
    //
    if (contains(params, "Materials")) {
      if (is<size_type>(params, "Materials")) {
        this->materials_identifiers.push_back(
            get<size_type>(params, "Materials"));
      } else {
        for (const auto& mid :
             get<std::vector<Parameter>>(params, "Materials")) {
          this->materials_identifiers.push_back(get<size_type>(mid));
        }
      }
    } else {
      this->materials_identifiers = p.getMaterialIdentifiers();
    }
    if (!contains(params, "Results")) {
      raise(
          "ParaviewExportIntegrationPointResultsAtNodes::"
          "ParaviewExportIntegrationPointResultsAtNodes: "
          "no results to export declared");
    }
    // 
    auto add_result = [this, &p](const std::string& rn) {
      auto r = IntegrationPointResult{};
      r.name = rn;
      this->getResultDescription(r, p);
      // creating the finite element space that will support the result grid
      // function
      auto fes = p.getFiniteElementSpace();
      if constexpr (parallel) {
        r.fespace = std::make_unique<FiniteElementSpace<parallel>>(
            fes.GetParMesh(), fes.FEColl(), r.number_of_components,
            fes.GetOrdering());
      } else {
        r.fespace = std::make_unique<FiniteElementSpace<parallel>>(
            fes.GetMesh(), fes.FEColl(), r.number_of_components,
            fes.GetOrdering());
      }
      // grid function
      r.f = std::make_unique<GridFunction<parallel>>(r.fespace.get());
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
      this->updateResultGridFunction(r, p);
    }
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  void ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      updateResultGridFunction(
          IntegrationPointResult& r,
          NonLinearEvolutionProblemImplementation<parallel>& p) {
    if (r.number_of_components == 1u) {
      auto c = ScalarResultCoefficient(r, p, this->materials_identifiers);
      r.f->ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    } else {
      auto c =
          MultiComponentsResultCoefficient(r, p, this->materials_identifiers);
      r.f->ProjectDiscCoefficient(c, mfem::GridFunction::ARITHMETIC);
    }
  }  // end of updateResultGridFunction

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
  ParaviewExportIntegrationPointResultsAtNodes<
      parallel>::~ParaviewExportIntegrationPointResultsAtNodes() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX */
