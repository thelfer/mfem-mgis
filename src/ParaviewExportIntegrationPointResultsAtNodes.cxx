/*!
 * \file   src/ParaviewExportIntegrationPointResultsAtNodes.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/05/2025
 */

#include <array>
#include <utility>
#include <optional>
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Profiler.hxx"
#ifdef MGIS_FUNCTION_SUPPORT
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/MaterialQuantityProviderSearch.hxx"
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
            if (f.getPartialQuadratureSpace().isDefinedOnABoundary()) {
              raise("functions defined on boundaries are not supported");
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
      attributes::Throwing,
      MaterialIntegrationPointResultBase& r,
      const NonLinearEvolutionProblemImplementationBase& p) {
    using namespace mgis::behaviour;
    using Status = MaterialQuantityProviderSearchResult::Status;
    using Category = MaterialIntegrationPointResultBase::Category;
    if (this->materials_identifiers.empty()) {
      raise("getResultDescription: empty list of material identifiers");
    }
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    auto c = Category{};
    auto s = size_type{};
    auto bis = std::vector<const AbstractBehaviourIntegrator*>{};
    bis.reserve(this->materials_identifiers.size());
    for (const auto& mid : this->materials_identifiers) {
      const auto l = LocationIdentifier{
          .material_identifier = MaterialIdentifier{.id = mid},
          .boundary_identifier = {}};
      const auto results =
          std::array<std::pair<Category, MaterialQuantityProviderSearchResult>,
                     3u>{
              std::pair{Category::GRADIENTS,
                        hasGradientProvider(ctx, p, l, r.name) | or_raise},
              std::pair{
                  Category::THERMODYNAMIC_FORCES,
                  hasThermodynamicForceProvider(ctx, p, l, r.name) | or_raise},
              std::pair{Category::INTERNAL_STATE_VARIABLES,
                        hasInternalStateVariableProvider(ctx, p, l, r.name) |
                            or_raise}};
      auto provider = std::optional<
          std::pair<Category, const AbstractBehaviourIntegrator*>>{};
      for (const auto& [c2, result] : results) {
        if (result.status == Status::NO_PROVIDER) {
          continue;
        }
        if ((result.status == Status::MULTIPLE_PROVIDERS) ||
            (provider.has_value())) {
          raise(
              "getResultDescription: multiple behaviour integrators "
              "provide the result '" +
              r.name + "' on material '" + std::to_string(mid) + "'");
        }
        provider = std::pair{c2, &(*(result.behaviour_integrator))};
      }
      if (!provider.has_value()) {
        raise(
            "getResultDescription: no behaviour integrator provides the "
            "result '" +
            r.name + "' on material '" + std::to_string(mid) + "'");
      }
      const auto [c2, bi] = *provider;
      const auto om = bi->getMaterial(ctx);
      if (isInvalid(om)) {
        raise(ctx.getErrorMessage());
      }
      const auto& m = *om;
      const auto& h = m.b.hypothesis;
      const auto& variables = [&m, c2]() -> const std::vector<Variable>& {
        if (c2 == Category::GRADIENTS) {
          return m.b.gradients;
        } else if (c2 == Category::THERMODYNAMIC_FORCES) {
          return m.b.thermodynamic_forces;
        }
        return m.b.isvs;
      }();
      const auto s2 = getVariableSize(getVariable(variables, r.name), h);
      if (bis.empty()) {
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
      bis.push_back(bi);
    }
    r.category = c;
    r.number_of_components = s;
    r.behaviour_integrators = std::move(bis);
  }  // end of getResultDescription

  std::optional<std::vector<ImmutablePartialQuadratureFunctionView>>
  ParaviewExportIntegrationPointResultsAtNodesBase::
      getPartialQuadratureFunctionViews(
          Context& ctx,
          const MaterialIntegrationPointResultBase& r,
          const TimeStepStage s) noexcept {
    using Category = MaterialIntegrationPointResultBase::Category;
    auto fcts = std::vector<ImmutablePartialQuadratureFunctionView>{};
    fcts.reserve(r.behaviour_integrators.size());
    for (const auto& bi : r.behaviour_integrators) {
      const auto om = bi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      const auto& m = *om;
      if (r.category == Category::GRADIENTS) {
        const auto og = getGradient(ctx, m, r.name, s);
        if (isInvalid(og)) {
          return {};
        }
        fcts.push_back(std::move(*og));
      } else if (r.category == Category::THERMODYNAMIC_FORCES) {
        const auto oth = getThermodynamicForce(ctx, m, r.name, s);
        if (isInvalid(oth)) {
          return {};
        }
        fcts.push_back(std::move(*oth));
      } else {
        const auto ov = getInternalStateVariable(ctx, m, r.name, s);
        if (isInvalid(ov)) {
          return {};
        }
        fcts.push_back(std::move(*ov));
      }
    }
    return fcts;
  }  // end of getPartialQuadratureFunctionViews

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(Context& ctx,
                                                   NonLinearEvolutionProblem& p,
                                                   const Parameters& params) {
    CatchTimeSection(ctx, "ParaviewExportResults::Constructor");
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          ctx, i, params);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          ctx, i, params);
    }
  }

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(
          Context& ctx,
          NonLinearEvolutionProblem& p,
          const ExportedFunctionsDescription& efcts,
          const std::string& d) {
    CatchTimeSection(ctx, "ParaviewExportResults::Constructor");
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          ctx, i, efcts, d);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          ctx, i, efcts, d);
    }
  }

  ParaviewExportIntegrationPointResultsAtNodes::
      ParaviewExportIntegrationPointResultsAtNodes(
          Context& ctx,
          NonLinearEvolutionProblem& p,
          const std::vector<ExportedFunctionsDescription>& ds,
          const std::string& d) {
    CatchTimeSection(ctx, "ParaviewExportResults::Constructor");
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
          ctx, i, ds, d);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto& i = p.getImplementation<false>();
      this->implementations.emplace<
          ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
          ctx, i, ds, d);
    }
  }

  bool ParaviewExportIntegrationPointResultsAtNodes::execute(
      Context& ctx,
      NonLinearEvolutionProblem& p,
      const real t,
      const real dt) noexcept {
    CatchTimeSection(ctx, "ParaviewExportResults::Execute");
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& i = p.getImplementation<true>();
      return std::
          get<ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
                 this->implementations)
              .execute(ctx, i, t, dt);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    auto& i = p.getImplementation<false>();
    return std::
        get<ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>(
               this->implementations)
            .execute(ctx, i, t, dt);
  }

  bool
  ParaviewExportIntegrationPointResultsAtNodes::executeInitialPostProcessing(
      Context& ctx, NonLinearEvolutionProblem& p, const real t) noexcept {
    CatchTimeSection(ctx,
                     "ParaviewExportIntegrationPointResultsAtNodes::"
                     "ExecuteInitialPostProcessing");
    const auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      using Implementation =
          ParaviewExportIntegrationPointResultsAtNodesImplementation<true>;
      auto& i = p.getImplementation<true>();
      return std::get<Implementation>(this->implementations)
          .executeInitialPostProcessing(ctx, i, t);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    using Implementation =
        ParaviewExportIntegrationPointResultsAtNodesImplementation<false>;
    auto& i = p.getImplementation<false>();
    return std::get<Implementation>(this->implementations)
        .executeInitialPostProcessing(ctx, i, t);
  }  // end of executeInitialPostProcessing

  ParaviewExportIntegrationPointResultsAtNodes::
      ~ParaviewExportIntegrationPointResultsAtNodes() = default;

#ifdef MGIS_FUNCTION_SUPPORT

  ParaviewExportIntegrationPointResultsAtNodesBase::ExportedFunctionsDescription
  makeExportedFunctionsDescription(std::string_view n,
                                   const PartialQuadratureFunctionsSet& f) {
    auto efcts = std::vector<ImmutablePartialQuadratureFunctionView>{};
    const auto& fcts = f.getFunctions();
    efcts.reserve(fcts.size());
    for (const auto& fptr : fcts) {
      efcts.push_back(*fptr);
    }
    return {.name = std::string{n}, .functions = std::move(efcts)};
  }

  std::vector<ParaviewExportIntegrationPointResultsAtNodesBase::
                  ExportedFunctionsDescription>
  makeExportedFunctionsDescriptions(
      const std::map<std::string, const PartialQuadratureFunctionsSet&>& fcts) {
    auto efcts = std::vector<ParaviewExportIntegrationPointResultsAtNodesBase::
                                 ExportedFunctionsDescription>{};
    for (const auto& [n, f] : fcts) {
      efcts.push_back(makeExportedFunctionsDescription(n, f));
    }
    return efcts;
  }

  PartialQuadratureFunctionsSet buildPartialQuadratureFunctionsSet(
      const NonLinearEvolutionProblemImplementationBase& p,
      const std::vector<size_type>& mids,
      const size_type nc) {
    auto qspaces = std::vector<std::shared_ptr<const PartialQuadratureSpace>>{};
    qspaces.reserve(mids.size());
    for (const auto m : mids) {
      qspaces.push_back(p.getMaterial(m).getPartialQuadratureSpacePointer());
    }
    return PartialQuadratureFunctionsSet(qspaces, nc);
  }  // end of buildPartialQuadratureFunctionsSet

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis