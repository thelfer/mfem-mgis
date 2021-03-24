/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export the results to paraview
   */
  template <bool parallel>
  struct ParaviewExportResults final
      : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportResults(NonLinearEvolutionProblemImplementation<parallel>&,
                          const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportResults() override;

   private:
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! exported grid function
    mfem_mgis::GridFunction<parallel> result;
  };  // end of struct ParaviewExportResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX */
