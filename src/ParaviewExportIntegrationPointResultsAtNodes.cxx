/*!
 * \file   src/ParaviewExportIntegrationPointResultsAtNodes.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/05/2025
 */

#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"

namespace mfem_mgis {

  ParaviewExportIntegrationPointResultsAtNodesBase::
      ParaviewExportIntegrationPointResultsAtNodesBase(const std::string& d)
      : exporter(d), cycle(0) {
    this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);
  }

  void
  ParaviewExportIntegrationPointResultsAtNodesBase::extractMaterialIdentifiers(
      const std::vector<ExportedFunctionsDescription>& ds) {
    auto get_materials_identifiers_and_check =
        [](const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) {
          auto mids = std::vector<size_type>{};
          if (fcts.empty()) {
            raise("no function given");
          }
          const auto nc = fcts.at(0).getNumberOfComponents();
          for (const auto& f : fcts) {
            if (f.getNumberOfComponents() != nc) {
              raise("inconsistent number of components");
            }
            const auto mid = f.getPartialQuadratureSpace().getId();
            if (std::find(mids.begin(), mids.end(), mid) != mids.end()) {
              raise("multiple function defined on material '" +
                    std::to_string(mid) + "'");
            }
            mids.push_back(mid);
          }
          return mids;
        };
    //
    if (ds.empty()) {
      raise("no functions given");
    }
    this->materials_identifiers =
        get_materials_identifiers_and_check(ds.at(0).functions);
    for (const auto& d : ds) {
      const auto mids = get_materials_identifiers_and_check(d.functions);
      if (mids.size() != this->materials_identifiers.size()) {
        raise("inconsistent material definitions");
      }
      if (!std::equal(this->materials_identifiers.begin(),
                      this->materials_identifiers.end(), mids.begin())) {
        raise("inconsistent material definitions");
      }
    }
  }

  void ParaviewExportIntegrationPointResultsAtNodesBase::getResultDescription(
      MaterialIntegrationPointResultBase& r,
      const NonLinearEvolutionProblemImplementationBase& p) {
    using namespace mgis::behaviour;
    auto c = MaterialIntegrationPointResultBase::Category{};
    auto s = size_type{};
    auto first = true;
    if (this->materials_identifiers.empty()) {
      raise("getResultDescription: empty list of material identifiers");
    }
    for (const auto& mid : this->materials_identifiers) {
      const auto [c2, s2] =
          [&p, &mid,
           &r]() -> std::tuple<MaterialIntegrationPointResultBase::Category,
                               size_type> {
        const auto& m = p.getMaterial(mid);
        const auto& h = m.b.hypothesis;
        if (contains(m.b.gradients, r.name)) {
          return {MaterialIntegrationPointResultBase::GRADIENTS,
                  getVariableSize(getVariable(m.b.gradients, r.name), h)};
        } else if (contains(m.b.thermodynamic_forces, r.name)) {
          return {MaterialIntegrationPointResultBase::THERMODYNAMIC_FORCES,
                  getVariableSize(getVariable(m.b.thermodynamic_forces, r.name),
                                  h)};
        } else if (!contains(m.b.isvs, r.name)) {
          raise("getResultDescription: no result '" + std::string(r.name) +
                "' found for material '" + std::to_string(mid) + "'");
        }
        return {MaterialIntegrationPointResultBase::INTERNAL_STATE_VARIABLES,
                getVariableSize(getVariable(m.b.isvs, r.name), h)};
      }();
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

  std::vector<ImmutablePartialQuadratureFunctionView>
  ParaviewExportIntegrationPointResultsAtNodesBase::
      getPartialQuadratureFunctionViews(
          const NonLinearEvolutionProblemImplementationBase& p,
          const MaterialIntegrationPointResultBase& r) {
    auto fcts = std::vector<ImmutablePartialQuadratureFunctionView>{};
    // creating the partial functions per materials
    for (const auto& mid : this->materials_identifiers) {
      const auto& m = p.getMaterial(mid);
      if (r.category == MaterialIntegrationPointResultBase::GRADIENTS) {
        fcts.push_back(getGradient(m, r.name));
      } else if (r.category ==
                 MaterialIntegrationPointResultBase::THERMODYNAMIC_FORCES) {
        fcts.push_back(getThermodynamicForce(m, r.name));
      } else {
        fcts.push_back(getInternalStateVariable(m, r.name));
      }
    }
    return fcts;
  }  // end of getPartialQuadratureFunctionViews

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(NonLinearEvolutionProblem& p,
                                                   const Parameters& params) {
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          i, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          i, params);
    }
  }

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblem& p,
          const ExportedFunctionsDescription& efcts,
          const std::string& d) {
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          i, efcts, d);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          i, efcts, d);
    }
  }

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblem& p,
          const std::vector<ExportedFunctionsDescription>& ds,
          const std::string& d) {
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          i, ds, d);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          i, ds, d);
    }
  }

  void ParaviewExportIntegrationPointResultsAtNodes::execute(
      NonLinearEvolutionProblem& p, const real t, const real dt) {
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      std::get<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          this->implementations)
          .execute(i, t, dt);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      std::get<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          this->implementations)
          .execute(i, t, dt);
    }
  }

  ParaviewExportIntegrationPointResultsAtNodes::
      ~ParaviewExportIntegrationPointResultsAtNodes() = default;

}  // end of namespace mfem_mgis
