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
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const Parameters& params)
      : ParaviewExportIntegrationPointResultsAtNodesBase(
            get<std::string>(params, "OutputFileName")) {
    checkParameters(params, {"OutputFileName", "Materials", "Results"});
    // if Materials exists, use it, otherwise, take all materials
    this->materials_identifiers = getMaterialsIdentifiers(p, params);
    this->createSubMesh(p);
    this->exporter.SetMesh(this->submesh.get());
    //
    if (!contains(params, "Results")) {
      raise(
          "ParaviewExportIntegrationPointResultsAtNodesImplementation::"
          "ParaviewExportIntegrationPointResultsAtNodesImplementation: "
          "no results to export declared");
    }
    //
    auto add_result = [this, &p](const std::string& rn) {
      auto r = MaterialIntegrationPointResult{};
      r.name = rn;
      this->getResultDescription(r, p);
      std::tie(r.fespace, r.f) = makeGridFunction<parallel>(
          this->getPartialQuadratureFunctionViews(p, r), this->submesh);
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
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      ParaviewExportIntegrationPointResultsAtNodesImplementation(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const ExportedFunctionsDescription& d,
          const std::string& n)
      : ParaviewExportIntegrationPointResultsAtNodesImplementation(
            p, std::vector<ExportedFunctionsDescription>(1, d), n) {
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      ParaviewExportIntegrationPointResultsAtNodesImplementation(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const std::vector<ExportedFunctionsDescription>& ds,
          const std::string& n)
      : ParaviewExportIntegrationPointResultsAtNodesBase(n) {
    this->extractMaterialIdentifiers(ds);
    this->createSubMesh(p);
    this->exporter.SetMesh(this->submesh.get());
    //
    for (const auto& d : ds) {
      auto fcts = std::make_unique<ExportedFunctions>();
      fcts->name = d.name;
      fcts->functions = d.functions;
      std::tie(fcts->grid_function_fespace, fcts->grid_function) =
          makeGridFunction<parallel>(fcts->functions, this->submesh);
      // registring
      this->exporter.RegisterField(fcts->name, fcts->grid_function.get());
      this->exported_functions.push_back(std::move(fcts));
    }
  }  // end of ParaviewExportIntegrationPointResultsAtNodesImplementation

  template <bool parallel>
  void ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::
      createSubMesh(NonLinearEvolutionProblemImplementation<parallel>& p) {
    /** Create Submesh using the material identifiers */
    mfem::Array<int> mat_attributes;
    for (const auto& mids : this->materials_identifiers) {
      mat_attributes.Append(mids);
    }
    this->submesh = std::make_shared<SubMesh<parallel>>(SubMesh<parallel>(
        SubMesh<parallel>::CreateFromDomain(p.getMesh(), mat_attributes)));
  }  // end of createSubMesh

  template <bool parallel>
  void
  ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);
    // updating grid functions
    if (!this->results.empty()) {
      for (auto& r : this->results) {
        updateGridFunction<parallel>(
            *(r.f), this->getPartialQuadratureFunctionViews(p, r),
            this->submesh);
      }
    } else {
      for (const auto& fcts : this->exported_functions) {
        updateGridFunction<parallel>(*(fcts->grid_function), fcts->functions,
                                     this->submesh);
      }
    }
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodesImplementation<
      parallel>::~ParaviewExportIntegrationPointResultsAtNodesImplementation() =
      default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_IXX */
