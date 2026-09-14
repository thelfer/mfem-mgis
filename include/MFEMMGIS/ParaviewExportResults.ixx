/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX

#include "mfem.hpp"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include <MFEMMGIS/Profiler.hxx>

namespace mfem_mgis {

  template <typename Mesh>
  void print_mesh_information(Mesh* mesh) {
    using Profiler::Utils::Message;
    using Profiler::Utils::sum;

    // get the number of vertices
    int64_t numbers_of_vertices_local = mesh->GetNV();
    int64_t numbers_of_vertices = sum(numbers_of_vertices_local);

    // get the number of elements
    int64_t numbers_of_elements_local = mesh->GetNE();
    int64_t numbers_of_elements = sum(numbers_of_elements_local);

    Message("INFO: number of vertices -> ", numbers_of_vertices);
    Message("INFO: number of elements -> ", numbers_of_elements);
  }

  template <bool parallel>
  ParaviewExportResults<parallel>::ParaviewExportResults(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& pb,
      const Parameters& params)
      : exporter(get<std::string>(throwing, params, "OutputFileName")),
        result(&pb.getFiniteElementSpace()),
        cycle(0) {
    //
    CatchTimeSection(ctx, "ParaviewExportResults::Constructor");
    //
    auto or_raise = ctx.getThrowingFailureHandler();
    //
    auto& u1 = pb.getUnknowns(ets);
    this->result.MakeTRef(&pb.getFiniteElementSpace(), u1, 0);

    /** Default, the mesh is the entire mesh */
    Mesh<parallel>& pmesh = pb.getMesh();

    bool contains_brd =
        contains(params, "Boundary") || contains(params, "Boundaries");
    bool contains_mat =
        contains(params, "Material") || contains(params, "Materials");

    if (contains_brd && contains_mat) {
      raise(
          "You can not define both materials and boundaries in a single "
          "ParaviewExportResults post processing");
    }
    if (contains_mat) {
      if (contains(params, "Material") && contains(params, "Materials")) {
        raise(
            "You can not use both 'Material' and 'Materials' parameters in a "
            "single ParaviewExportResults post processing");
      }
    }
    if (contains_brd) {
      if (contains(params, "Boundary") && contains(params, "Boundaries")) {
        raise(
            "You can not use both 'Boundary' and 'Boundaries' parameters in a "
            "single ParaviewExportResults post processing");
      }
    }
    //
    if (contains_brd || contains_mat) {
      const auto l = [&contains_mat] {
        if (contains_mat) {
          return MeshDiscretization::Location::ON_MATERIALS;
        }
        return MeshDiscretization::Location::ON_BOUNDARIES;
      }();
      const auto ids = [&contains_mat, &contains_brd, &params] {
        if (contains_mat) {
          return get(throwing, params,
                     contains(params, "Material") ? "Material" : "Materials");
        }
        return get(throwing, params,
                   contains(params, "Boundary") ? "Boundary" : "Boundaries");
      }();
      const auto on_all_materials = [&contains_mat, &pb, &ids] {
        if (contains_mat) {
          const auto& fed = pb.getFiniteElementDiscretization();
          const auto mids = getMaterialsIdentifiers(throwing, fed, ids);
          return static_cast<size_type>(mids.size()) ==
                 getMaterialsAttributes(fed).Size();
        }
        return false;
      }();
      if (!on_all_materials) {
        const auto& fed = pb.getFiniteElementDiscretization();
        this->submesh =
            fed.template getMutableSubMeshPointer<parallel>(ctx, ids, l) |
            or_raise;
        auto fespaces_manager =
            pb.getFiniteElementDiscretization().getFiniteElementSpacesManager();
        const auto nc =
            fed.template getFiniteElementSpace<parallel>().GetVDim();
        /** create the underlying finite element space */
        this->fes_sm =
            fespaces_manager.template getFiniteElementSpace<parallel>(
                ctx, {.location = l,
                      .identifiers = ids,
                      .number_of_components = nc}) |
            or_raise;
        /** init the grid function corresponding to the sub mesh */
        this->result_sm =
            std::make_shared<mfem_mgis::GridFunction<parallel>>(fes_sm.get());
      }
    }
    // setting the exporter
    this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);
    if (this->submesh.get() != nullptr) {
      /** Update exporter */
      this->exporter.SetMesh(this->submesh.get());
      if (contains(params, "Verbosity")) {
        if (get<int>(throwing, params, "Verbosity") >= 1) {
          Profiler::Utils::Message(
              "Submesh information [for domain attributes]");
          print_mesh_information(this->submesh.get());
        }
      }
      if (contains(params, "OutputFieldName")) {
        this->exporter.RegisterField(
            get<std::string>(throwing, params, "OutputFieldName"),
            this->result_sm.get());
      } else {
        this->exporter.RegisterField("u", this->result_sm.get());
      }
    } else { /** Not a sub mesh */
      exporter.SetMesh(&pmesh);
      if (contains(params, "OutputFieldName")) {
        this->exporter.RegisterField(
            get<std::string>(throwing, params, "OutputFieldName"),
            &(this->result));
      } else {
        this->exporter.RegisterField("u", &(this->result));
      }
    }
  }  // end of ParaviewExportResults

  template <bool parallel>
  void ParaviewExportResults<parallel>::execute(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>&,
      const real t,
      const real dt) {
    CatchTimeSection(ctx, "ParaviewExportResults::Execute");
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);
    // SetFromTrueVector needed here in MFEM for at least two rationales:
    //    - it applies prolongation matrix (Non-Conforming mesh, BCs, AMR ...)
    //      to set the values of some unkwown dofs deduced from known dofs
    //    - exchange data between processes in order to retrieve information
    //      needed to perform the previous prolongation step
    if (submesh.get() != nullptr) {
      /** Transfer data from global mesh to submesh */
      this->result.SetFromTrueVector();
      this->submesh.get()->Transfer(this->result, this->result_sm.get()[0]);
      this->exporter.Save();
    } else {
      this->result.SetFromTrueVector();
      this->exporter.Save();
    }
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
