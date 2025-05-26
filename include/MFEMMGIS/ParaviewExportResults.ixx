/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX

#include "mfem.hpp"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include <MFEMMGIS/Profiler.hxx>

namespace mfem_mgis {

  template<typename Mesh>
    void print_mesh_information(Mesh* mesh)
    {

      using Profiler::Utils::sum;
      using Profiler::Utils::Message;

      //get the number of vertices
      int64_t numbers_of_vertices_local = mesh->GetNV();
      int64_t  numbers_of_vertices = sum(numbers_of_vertices_local);

      //get the number of elements
      int64_t numbers_of_elements_local = mesh->GetNE();
      int64_t numbers_of_elements = sum(numbers_of_elements_local);

      Message("INFO: number of vertices -> ", numbers_of_vertices);
      Message("INFO: number of elements -> ", numbers_of_elements);
    }

  template <bool parallel>
    ParaviewExportResults<parallel>::ParaviewExportResults(
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const Parameters& params)
    : exporter(get<std::string>(params, "OutputFileName")),
    result(&p.getFiniteElementSpace()),
    cycle(0) {
      CatchTimeSection("ParaviewExportResults::Constructor");
      auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
      this->result.MakeTRef(&p.getFiniteElementSpace(), u1, 0);

      /** Default, the mesh is the entire mesh */
      Mesh<parallel>& pmesh = p.getMesh();
      
      bool contains_brd = contains(params, "Boundary") || contains(params, "Boundaries");
      bool contains_mat = contains(params, "Material") || contains(params, "Materials");

      if(contains_brd && contains_mat)
      {
        Profiler::Utils::Message("You can not define 'Material' and 'Boundary' in a single ParaviewExportResults post processing");
        std::exit(0);
      }

      if(contains_mat){ /** Materials and Sub mesh */
        /** "false" means that we double check if params include Material or Materials */
        auto materials_ids = getMaterialsIdentifiers(p, params, false);
        mfem::Array<int> mat_attributes;

        /** Create Submesh using the material identifiers */
        for(const auto& mids : materials_ids){
          mat_attributes.Append(mids);
        }
        
        this->submesh = std::make_shared<SubMesh<parallel>>(
        SubMesh<parallel>(
        SubMesh<parallel>::CreateFromDomain(pmesh, mat_attributes)
        ));

        /** Create the corresponding Grid Function */
        auto& FED = p.getFiniteElementDiscretization();
        const auto& fec_subdomain = FED.getFiniteElementCollection();

        this->fes_sm = std::make_shared<FiniteElementSpace<parallel>>(FiniteElementSpace<parallel>( 
              this->submesh.get(), 
              &fec_subdomain, 
              FED.template getFiniteElementSpace<parallel>().GetVDim()));

        /** init the grid function corresponding to the sub mesh */
        this->result_sm = std::make_shared<mfem_mgis::GridFunction<parallel>>(
        mfem_mgis::GridFunction<parallel>(fes_sm.get()));

        /** Update exporter */
        this->exporter.SetMesh(this->submesh.get());
        this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);

        if (contains(params, "Verbosity")) {
          if(get<int>(params, "Verbosity") >= 1) {
            Profiler::Utils::Message("Submesh information [for domain attributes]");
            print_mesh_information(this->submesh.get());
          }
        }

        if (contains(params, "OutputFieldName")) {
          this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
              this->result_sm.get());
        } else {
          this->exporter.RegisterField("u", this->result_sm.get());
        }
      } /** No Domain Attributes */
      else if(contains_brd){ /** Boundaries and Sub mesh */
      
        /** "false" means that we double check if params include Material or Materials */
        auto bdrAttributes = getBoundariesIdentifiers(p , params, false);

        /** Get the list of boundary attributes */
        mfem::Array<int> bdr_attributes;
        
        /** Create Submesh using the domain attributes */
        for(const auto& dattr : bdrAttributes){
          bdr_attributes.Append(dattr);
        }

        /** Use the list of boundary attributes to define the sub mesh */
        this->submesh = std::make_shared<SubMesh<parallel>>(
        SubMesh<parallel>(
        SubMesh<parallel>::CreateFromBoundary(pmesh, bdr_attributes)));

        /** Create the corresponding Grid Function */
        auto& FED = p.getFiniteElementDiscretization();
        const auto& fec_subdomain = FED.getFiniteElementCollection();

        this->fes_sm = std::make_shared<FiniteElementSpace<parallel>>(FiniteElementSpace<parallel>( 
              this->submesh.get(), 
              &fec_subdomain, 
              pmesh.Dimension()));

        /** init the grid function corresponding to the sub mesh */
        this->result_sm = std::make_shared<mfem_mgis::GridFunction<parallel>>(
        mfem_mgis::GridFunction<parallel>(fes_sm.get()));

        /** Update exporter */
        this->exporter.SetMesh(this->submesh.get());
        this->exporter.SetDataFormat(mfem::VTKFormat::BINARY);

        if (contains(params, "Verbosity")) {
          if(get<int>(params, "Verbosity") >= 1) {
            Profiler::Utils::Message("Submesh information [for boundary attributes]");
            print_mesh_information(this->submesh.get());
          }
        }

        if (contains(params, "OutputFieldName")) {
          this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
              this->result_sm.get());
        } else {
          this->exporter.RegisterField("u", this->result_sm.get());
        }
      } /** No BoundaryAttributes */
      else { /** Not a sub mesh */
        exporter.SetMesh(&pmesh);
        if (contains(params, "OutputFieldName")) {
          this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
              &(this->result));
        } else {
          this->exporter.RegisterField("u", &(this->result));
        }
      }
    }  // end of ParaviewExportResults

  template <bool parallel>
    void ParaviewExportResults<parallel>::execute(
        NonLinearEvolutionProblemImplementation<parallel>&,
        const real t,
        const real dt) {
      CatchTimeSection("ParaviewExportResults::Execute");
      this->exporter.SetCycle(this->cycle);
      this->exporter.SetTime(t + dt);
      // SetFromTrueVector needed here in MFEM for at least two rationales:
      //    - it applies prolongation matrix (Non-Conforming mesh, BCs, AMR ...)
      //      to set the values of some unkwown dofs deduced from known dofs
      //    - exchange data between processes in order to retrieve information
      //      needed to perform the previous prolongation step
      if(submesh != nullptr)
      {
        /** Transfer data from gloabl mesh to submesh */
        this->result.SetFromTrueVector();
        this->submesh.get()->Transfer(this->result, this->result_sm.get()[0]);
        this->exporter.Save();
      }
      else
      {
        this->result.SetFromTrueVector();
        this->exporter.Save();
      }
      ++(this->cycle);
    }  // end of execute

  template <bool parallel>
    ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
