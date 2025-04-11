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

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportResults<parallel>::ParaviewExportResults(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        result(&p.getFiniteElementSpace()),
        cycle(0) {
    auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    this->result.MakeTRef(&p.getFiniteElementSpace(), u1, 0);
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->result));
    } else if(contains(params, "DomainAttributes")){
          auto domainAttributes = get<std::vector<int>>(params, "DomainAttributes");
          for(const auto& dattr : domainAttributes){
          	this->domain_attributes.Append(dattr);
          }
    } else {
      this->exporter.RegisterField("u", &(this->result));
    }
  }  // end of ParaviewExportResults

  template <bool parallel>
  void ParaviewExportResults<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);
    // SetFromTrueVector needed here in MFEM for at least two rationales:
    //    - it applies prolongation matrix (Non-Conforming mesh, BCs, AMR ...)
    //      to set the values of some unkwown dofs deduced from known dofs
    //    - exchange data between processes in order to retrieve information
    //      needed to perform the previous prolongation step
    if(this->domain_attributes.Size() != 0)
    {
    std::cout << "with domain attribute " << std::endl;
    	/** Version with copies */
    	Mesh<parallel>& pmesh = p.getMesh();
	mfem::Array<int> DomAtt = mfem::Array<int>(this->domain_attributes.Size());
	int it = 0;
	for(const auto& AttrType : this->domain_attributes) {
	   assert(AttrType <= pmesh.attributes.Max());
	   DomAtt[it++] = AttrType;
	}
	
	/** Get Sub mesh */
	SubMesh<parallel> PSubMesh = SubMesh<parallel>::CreateFromDomain(pmesh, DomAtt);
	
	/** Project data */
	auto& FED = p.getFiniteElementDiscretization();
	const auto& fec_subdomain = FED.getFiniteElementCollection();
	FiniteElementSpace<parallel> FESpaceSubMesh( &PSubMesh, &fec_subdomain, pmesh.Dimension());
	mfem_mgis::GridFunction<parallel> resultOnSubMesh(&FESpaceSubMesh);
	PSubMesh.Transfer(this->result, resultOnSubMesh);
	
	/** Export */
	mfem::ParaViewDataCollection subMeshExporter("test", &PSubMesh);
	subMeshExporter.SetCycle(this->cycle);
	subMeshExporter.SetTime(t + dt);
	subMeshExporter.RegisterField("Displacement", &resultOnSubMesh);
	subMeshExporter.SetDataFormat(mfem::VTKFormat::BINARY);
	subMeshExporter.Save();
    }
    else
    {
        std::cout << "without domain attribute " << std::endl;
    	this->result.SetFromTrueVector();
    	this->exporter.Save();
    }
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
