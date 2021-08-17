/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "mfem/fem/gridfunc.hpp"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportResults<parallel>::ParaviewExportResults(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        displacement(&p.getFiniteElementSpace()),
	nb_materials(p.getFiniteElementSpace().GetMesh()->attributes.Size()),
	dim(p.getFiniteElementSpace().GetMesh()->Dimension()),
	sdim(dim*(dim+1)/2),
        cycle(0) {
    auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    this->displacement.MakeTRef(&p.getFiniteElementSpace(), u1, 0);
    this->displacement.SetFromTrueVector();
    // Set up the list of materials
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->displacement));
    } else {
      this->exporter.RegisterField("u", &(this->displacement));
    }

    auto fes  = &(p.getFiniteElementSpace());
    auto fec  = fes->FEColl();
    auto ordering  = fes->GetOrdering();
    // scalarfield defines a scalar field type that covers the same FE collection that
    // the one of the problem (original FE space can be a vector field)
    mfem_mgis::FiniteElementSpace<parallel>* scalarfield;
    if constexpr (parallel) {
      scalarfield = new mfem_mgis::FiniteElementSpace<parallel>(fes->GetParMesh(),fec,1,ordering);
    } else {
      scalarfield = new mfem_mgis::FiniteElementSpace<parallel>(fes->GetMesh(),fec,1,ordering);
    }
    // stress_c and strain_c will be used for temporary calculations for each material
    stress_c.resize(nb_materials);
    strain_c.resize(nb_materials);
    mgis_materials.resize(nb_materials);
    stress.resize(sdim*nb_materials);
    strain.resize(sdim*nb_materials);   
    // Loop on materials
    for (mfem_mgis::size_type sm = 0; sm != nb_materials; ++sm) {
      // pointer on the considered material
      mgis_materials[sm] = &p.getMaterial(sm+1);
      stress_c[sm]= new StressCoefficient(mgis_materials[sm]);
      strain_c[sm]= new StrainCoefficient(mgis_materials[sm]);

      // Loop over the different components
      int loc_index = 0;
      for (int si = 0; si < dim; si++) {
	for (int sj = si; sj < dim; sj++, loc_index++) {
	  std::string letters = "xyz";
	  int index = sdim*sm + loc_index;
	  // Build a name for stress diagnostic
	  std::string stressname = "S" + letters.substr(si,1) +
	    letters.substr(sj,1) + "_mat_" + std::to_string(sm+1);
	  // Build a name for strain diagnostic
	  std::string strainname = (si != sj ? "G" : "E") +
	    letters.substr(si,1) + letters.substr(sj,1) +
	    "_mat_" + std::to_string(sm+1);
	  // Build GridFunction that will store the diagnostic
	  std::cout << "index " << index << " max " <<  sdim*nb_materials << std::endl;

	  stress[index] = new mfem_mgis::GridFunction<parallel>(scalarfield);
	  strain[index] = new mfem_mgis::GridFunction<parallel>(scalarfield);
	  // Declare strain and stress associated with material mat
	  // and component (si,sj)
	  this->exporter.RegisterField(stressname,stress[index]);
	  this->exporter.RegisterField(strainname,strain[index]);
	}
      }
    }
  }  // end of ParaviewExportResults

  template <bool parallel>
  void ParaviewExportResults<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>&,
      const real t,
      const real dt) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);
    this->displacement.SetFromTrueVector();

    // Loop on materials
    for (mfem_mgis::size_type sm = 0; sm != nb_materials; ++sm) {
      // Loop over the different components
      int loc_index = 0;
      for (int si = 0; si < dim; si++) {
	for (int sj = si; sj < dim; sj++, loc_index++) {
	  int index = sdim*sm + loc_index;

	  // Set component to be computed in temporary objects stress_c[m] and strain_c[m]
	  stress_c[sm]->SetComponent(si, sj);
	  strain_c[sm]->SetComponent(si, sj);
	  // Update displacement within temporary objects stress_c[m] and strain_c[m]
	  stress_c[sm]->SetDisplacement(this->displacement);
	  strain_c[sm]->SetDisplacement(this->displacement);
//	  // Perform the projection on stress[index] and strain[index]
//	  stress[index]->ProjectDiscCoefficient(*(stress_c[sm]), mfem::GridFunction::ARITHMETIC);
//	  strain[index]->ProjectDiscCoefficient(*(strain_c[sm]), mfem::GridFunction::ARITHMETIC);
	}
      }
    }

    // Export all registered fields
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
