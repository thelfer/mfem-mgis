/*!
 * \file   include/MFEMMGIS/ParaviewExportFullFieldResults.hxx
 * \brief
 * \author Maxence Wangermez
 * \date   10/10/2023
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_HXX

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export the results to paraview
   */
  template <bool parallel>
  struct ParaviewExportFullFieldResults final : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportFullFieldResults(NonLinearEvolutionProblemImplementation<parallel>&,
                          const Parameters&);
    
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;

    // mfem_mgis::real* getCoordinates(const mfem::GridFunction&,
    //                     const bool,
    //                     const size_t,
    //                     const int,
    //                     const int);
    
    //! \brief destructor
    ~ParaviewExportFullFieldResults() override;

   private:
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief exported grid function
    mfem_mgis::GridFunction<parallel> result;
    //! \brief number of records
    size_type cycle;
    //! \brief compute the fullfield solution
    void computeFullFieldResult(real&, const real&, const real&);
  };  // end of struct ParaviewExportFullFieldResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportFullFieldResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_HXX */
