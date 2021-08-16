/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportResults<parallel>::ParaviewExportResults(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        displacement(&p.getFiniteElementSpace()),
	nb_materials(p.getFiniteElementSpace().GetMesh()->attributes.Size()),
	mgis_materials(nb_materials),
	dim(p.getFiniteElementSpace().GetMesh()->Dimension()),
	sdim(dim*(dim+1)/2),
        cycle(0),
	stress(sdim),
	strain(sdim) {
    auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    this->displacement.MakeTRef(&p.getFiniteElementSpace(), u1, 0);
    this->displacement.SetFromTrueVector();
    // Set up the list of materials
    for (mfem_mgis::size_type i = 0; i != nb_materials; ++i) {
      mgis_materials[i] = &p.getMaterial(i);
    }
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->displacement));
    } else {
      this->exporter.RegisterField("u", &(this->displacement));
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
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
