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
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName")),
        cycle(0) {
    checkParameters(params, {"OutputFileName", "Materials", "Results"});
    // if Materials exists, use it, otherwise, take all materials
    this->materials_identifiers = getMaterialsIdentifiers(p, params);
    /** Create Submesh using the material identifiers */
    mfem::Array<int> mat_attributes;
    for (const auto& mids : this->materials_identifiers) {
      mat_attributes.Append(mids);
    }
    this->submesh = std::make_shared<SubMesh<parallel>>(SubMesh<parallel>(
        SubMesh<parallel>::CreateFromDomain(p.getMesh(), mat_attributes)));
    this->exporter.SetMesh(this->submesh.get());
    this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);
    //
    if (!contains(params, "Results")) {
      raise(
          "ParaviewExportIntegrationPointResultsAtNodes::"
          "ParaviewExportIntegrationPointResultsAtNodes: "
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
  }  // end of ParaviewExportIntegrationPointResultsAtNodes

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const ExportedFunctionsDescription& d,
          const std::string& n)
      : ParaviewExportIntegrationPointResultsAtNodes(
            p, std::vector<ExportedFunctionsDescription>(1, d), n) {
  }  // end of ParaviewExportIntegrationPointResultsAtNodes

  template <bool parallel>
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::
      ParaviewExportIntegrationPointResultsAtNodes(
          NonLinearEvolutionProblemImplementation<parallel>& p,
          const std::vector<ExportedFunctionsDescription>& ds,
          const std::string& n)
      : exporter(n), cycle(0) {
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
    for (const auto& d: ds) {
      const auto mids = get_materials_identifiers_and_check(d.functions);
      if (mids.size() != this->materials_identifiers.size()) {
        raise("inconsistent material definitions");
      }
      if (!std::equal(this->materials_identifiers.begin(),
                      this->materials_identifiers.end(), mids.begin())) {
        raise("inconsistent material definitions");
      }
    }
    /** Create Submesh using the material identifiers */
    mfem::Array<int> mat_attributes;
    for (const auto& mids : this->materials_identifiers) {
      mat_attributes.Append(mids);
    }
    this->submesh = std::make_shared<SubMesh<parallel>>(SubMesh<parallel>(
        SubMesh<parallel>::CreateFromDomain(p.getMesh(), mat_attributes)));
    this->exporter.SetMesh(this->submesh.get());
    this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);
    //
    for (const auto& d : ds) {
      auto fcts = std::make_unique<ExportedFunctions>();
      fcts->name = d.name; fcts->functions = d.functions;
      std::tie(fcts->grid_function_fespace, fcts->grid_function) =
          makeGridFunction<parallel>(fcts->functions, this->submesh);
      // registring
      this->exporter.RegisterField(fcts->name, fcts->grid_function.get());
      this->exported_functions.push_back(std::move(fcts));
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
    if (!this->results.empty()) {
      for (auto& r : this->results) {
        updateGridFunction<parallel>(
            *(r.f), this->getPartialQuadratureFunctionViews(p, r),
            this->submesh);
      }
    } else {
      for (const auto& fcts : this->exported_functions){
        updateGridFunction<parallel>(*(fcts->grid_function), fcts->functions,
                                     this->submesh);
      }
    }
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  void
  ParaviewExportIntegrationPointResultsAtNodes<parallel>::getResultDescription(
      MaterialIntegrationPointResult& r,
      const NonLinearEvolutionProblemImplementation<parallel>& p) {
    using namespace mgis::behaviour;
    auto c = typename MaterialIntegrationPointResult::Category{};
    auto s = size_type{};
    auto first = true;
    if (this->materials_identifiers.empty()) {
      raise("getResultDescription: empty list of material identifiers");
    }
    for (const auto& mid : this->materials_identifiers) {
      const auto [c2, s2] = [
        &p, &mid, &r
      ]() -> std::tuple<typename MaterialIntegrationPointResult::Category, size_type> {
        const auto& m = p.getMaterial(mid);
        const auto& h = m.b.hypothesis;
        if (contains(m.b.gradients, r.name)) {
          return {MaterialIntegrationPointResult::GRADIENTS,
                  getVariableSize(getVariable(m.b.gradients, r.name), h)};
        } else if (contains(m.b.thermodynamic_forces, r.name)) {
          return {MaterialIntegrationPointResult::THERMODYNAMIC_FORCES,
                  getVariableSize(getVariable(m.b.thermodynamic_forces, r.name),
                                  h)};
        } else if (!contains(m.b.isvs, r.name)) {
          raise("getResultDescription: no result '" + std::string(r.name) +
                "' found for material '" + std::to_string(mid) + "'");
        }
        return {MaterialIntegrationPointResult::INTERNAL_STATE_VARIABLES,
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
          const MaterialIntegrationPointResult& r) {
    auto fcts = std::vector<ImmutablePartialQuadratureFunctionView>{};
    // creating the partial functions per materials
    for (const auto& mid : this->materials_identifiers) {
      const auto& m = p.getMaterial(mid);
      if (r.category == MaterialIntegrationPointResult::GRADIENTS) {
        fcts.push_back(getGradient(m, r.name));
      } else if (r.category == MaterialIntegrationPointResult::THERMODYNAMIC_FORCES) {
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
