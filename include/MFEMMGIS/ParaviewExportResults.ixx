/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportResults<parallel>::ParaviewExportResults(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params)
      : exporter(get<std::string>(params, "FileName"),
                 p.getFiniteElementSpace().GetMesh()),
        result(&p.getFiniteElementSpace()) {
    auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    this->result.MakeTRef(&p.getFiniteElementSpace(), u1, 0);
    this->result.SetFromTrueVector();
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->result));
    } else {
      this->exporter.RegisterField("u", &(this->result));
    }
  }  // end of ParaviewExportResults

  template <bool parallel>
  void ParaviewExportResults<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>&,
      const real t,
      const real dt) {
    this->exporter.SetCycle(0);
    this->exporter.SetTime(t + dt);
    this->exporter.Save();
  }  // end of execute

  template <bool parallel>
  ParaviewExportResults<parallel>::~ParaviewExportResults() = default;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_IXX */
