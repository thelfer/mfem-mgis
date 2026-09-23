/*!
 * \file   include/MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.ixx
 * \brief
 * \author Thomas Helfer
 * \date   18/08/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX

#include "mfem/mesh/submesh/submesh.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/submesh/psubmesh.hpp"
#endif /* MFEM_USE_MPI */

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      ParaviewExportIntegrationPointResultsAtNodesImplementation(
          Context& ctx,
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const Parameters& params)
      : ParaviewExportIntegrationPointResultsAtNodesBase(
            get<std::string>(throwing, params, "OutputFileName")),
        shallExecuteInitialPostProcessing(get_if<bool>(
            throwing, params, "ExecuteInitialPostProcessing", true)) {
    checkParameters(throwing, params,
                    {"OutputFileName", "Materials", "Results",
                     "ExecuteInitialPostProcessing"});
    // if Materials exists, use it, otherwise, take all materials
    this->materials_identifiers = getMaterialsIdentifiers(throwing, p, params);
    const auto all_mids = p.getAssignedMaterialsIdentifiers();
    if (this->materials_identifiers.size() == all_mids.size()) {
      this->exporter.SetMesh(&(p.getMesh()));
    } else {
      this->createSubMesh(ctx, p);
      this->exporter.SetMesh(this->submesh.get());
    }
    //
    if (!contains(params, "Results")) {
      raise(
          "ParaviewExportIntegrationPointResultsAtNodesImplementation::"
          "ParaviewExportIntegrationPointResultsAtNodesImplementation: "
          "no results to export declared");
    }
    //
    auto add_result = [this, &p, &ctx](const std::string& rn) {
      auto r = MaterialIntegrationPointResult{};
      r.name = rn;
      this->getResultDescription(throwing, r, p);
      if (this->submesh.get() == nullptr) {
        auto or2 = makeGridFunction<parallel>(
            ctx, this->getPartialQuadratureFunctionViews(throwing, r),
            p.getMesh());
        if (isInvalid(or2)) {
          raise(ctx.getErrorMessage());
        }
        r.f = std::move(or2);
      } else {
        auto or2 = makeGridFunction<parallel>(
            ctx, this->getPartialQuadratureFunctionViews(throwing, r),
            *(this->submesh));
        if (isInvalid(or2)) {
          raise(ctx.getErrorMessage());
        }
        r.f = std::move(or2);
      }
      // registring
      this->exporter.RegisterField(r.name, r.f.get());
      // saving
      this->results.push_back(std::move(r));
    };
    if (is<std::string>(throwing, params, "Results")) {
      add_result(get<std::string>(throwing, params, "Results"));
    } else {
      for (const auto& rn :
           get<std::vector<Parameter>>(throwing, params, "Results")) {
        add_result(get<std::string>(throwing, rn));
      }
    }
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      ParaviewExportIntegrationPointResultsAtNodesImplementation(
          Context& ctx,
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const ExportedFunctionsDescription& d,
          const std::string& n)
      : ParaviewExportIntegrationPointResultsAtNodesImplementation(
            ctx, p, std::vector<ExportedFunctionsDescription>(1, d), n) {
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      ParaviewExportIntegrationPointResultsAtNodesImplementation(
          Context& ctx,
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const std::vector<ExportedFunctionsDescription>& ds,
          const std::string& n)
      : ParaviewExportIntegrationPointResultsAtNodesBase(n),
        shallExecuteInitialPostProcessing(true) {
    this->extractMaterialIdentifiers(ds);
    const auto all_mids = p.getAssignedMaterialsIdentifiers();
    if (this->materials_identifiers.size() == all_mids.size()) {
      this->exporter.SetMesh(&(p.getMesh()));
    } else {
      this->createSubMesh(ctx, p);
      this->exporter.SetMesh(this->submesh.get());
    }
    for (const auto& d : ds) {
      auto fcts = std::make_unique<ExportedFunctions>();
      fcts->name = d.name;
      fcts->functions = d.functions;
      if (this->submesh.get() == nullptr) {
        auto ores =
            makeGridFunction<parallel>(ctx, fcts->functions, p.getMesh());
        if (isInvalid(ores)) {
          raise(ctx.getErrorMessage());
        }
        fcts->grid_function = std::move(ores);
      } else {
        auto ores =
            makeGridFunction<parallel>(ctx, fcts->functions, *(this->submesh));
        if (isInvalid(ores)) {
          raise(ctx.getErrorMessage());
        }
        fcts->grid_function = std::move(ores);
      }
      // registring
      this->exporter.RegisterField(fcts->name, fcts->grid_function.get());
      this->exported_functions.push_back(std::move(fcts));
    }
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  void ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      createSubMesh(Context& ctx,
                    NonLinearEvolutionProblemImplementation<parallel>& p) {
    auto or_raise = ctx.getThrowingFailureHandler();
    auto fed = p.getFiniteElementDiscretization();
    this->submesh = fed.template getMutableSubMeshPointer<parallel>(
                        ctx, Parameter::from(this->materials_identifiers),
                        MeshDiscretization::Location::ON_MATERIALS) |
                    or_raise;
  }  // end of createSubMesh

  template <bool parallel>
  bool ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      executeInitialPostProcessing(
          Context&,
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const real t) noexcept {
    if (this->shallExecuteInitialPostProcessing) {
      this->exportResults(p, t, bts);
    }
    return true;
  }  // end of executeInitialPostProcessing

  template <bool parallel>
  void
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::execute(
      Context&,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    this->exportResults(p, t + dt, ets);
  }  // end of execute

  template <bool parallel>
  void ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      exportResults(NonLinearEvolutionProblemImplementation<parallel>& p,
                    const real t,
                    const TimeStepStage s) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t);
    // updating grid functions
    if (!this->results.empty()) {
      for (auto& r : this->results) {
        if (this->submesh.get() == nullptr) {
          updateGridFunction<parallel>(
              *(r.f), this->getPartialQuadratureFunctionViews(throwing, r, s),
              p.getMesh());
        } else {
          updateGridFunction<parallel>(
              *(r.f), this->getPartialQuadratureFunctionViews(throwing, r, s),
              *(this->submesh));
        }
      }
    } else {
      for (const auto& fcts : this->exported_functions) {
        if (this->submesh.get() == nullptr) {
          updateGridFunction<parallel>(*(fcts->grid_function), fcts->functions,
                                       p.getMesh());
        } else {
          updateGridFunction<parallel>(*(fcts->grid_function), fcts->functions,
                                       *(this->submesh));
        }
      }
    }
    this->exporter.Save();
    ++(this->cycle);
  }  // end of exportResults

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<
      parallel>::~ParaviewExportIntegrationPointResultsAtNodesImplementation() =
      default;

#ifdef MGIS_FUNCTION_SUPPORT

  template <bool parallel>
  ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointPostProcessingsResultsAtNodes(
          Context& ctx,
          NonLinearEvolutionProblemImplementation<parallel>& p,
          std::string_view n,
          const std::vector<size_type> mids,
          const size_type nc,
          std::function<bool(Context&, PartialQuadratureFunction&)> f,
          std::string_view d)
      : functions(buildPartialQuadratureFunctionsSet(p, mids, nc)),
        update_function(f),
        exporter(ctx,
                 p,
                 makeExportedFunctionsDescription(n, this->functions),
                 std::string{d}) {
  }  // end of
     // ParaviewExportIntegrationPointPostProcessingsResultsAtNodes

  template <bool parallel>
  bool ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<parallel>::
      executeInitialPostProcessing(
          Context&,
          NonLinearEvolutionProblemImplementation<parallel>&,
          const real) noexcept {
    return true;
  }  // end of executeInitialPostProcessing

  template <bool parallel>
  void ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
      parallel>::execute(Context& ctx,
                         NonLinearEvolutionProblemImplementation<parallel>& p,
                         const real t,
                         const real dt) {
    Context local_ctx;
    if (!this->functions.update(local_ctx, this->update_function)) {
      raise(ctx.getErrorMessage());
    }
    this->exporter.execute(ctx, p, t, dt);
  }  // end of execute

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX */
