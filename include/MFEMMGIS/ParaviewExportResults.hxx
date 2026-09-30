/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.hxx
 * \brief  This file declares the `ParaviewExportResults` class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export the results to paraview
   */
  template <bool parallel>
  struct ParaviewExportResults final
      : public AbstractNonLinearEvolutionProblemPostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] pb: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportResults(mgis::Context& ctx,
                          NonLinearEvolutionProblemImplementation<parallel>& pb,
                          const Parameters& params);
    /*!
     * \brief execute the post-processing at the initial time
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t) noexcept override;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem, unused
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t,
        const real dt) noexcept override;
    //! \brief destructor
    ~ParaviewExportResults() override;

   private:
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief exported grid function
    mfem_mgis::GridFunction<parallel> result;
    //! \brief number of records
    size_type cycle;
    /*!
     * \brief boolean stating if the results shall be exported at the initial
     * time of the simulation
     */
    const bool shallExecuteInitialPostProcessing;
    //! \brief submesh defined when exporting data for domain or boundary
    //! attributes
    std::shared_ptr<mfem_mgis::SubMesh<parallel>> submesh;
    //! \brief fespace defined when exporting data for domain or boundary
    //! attributes
    std::shared_ptr<mfem_mgis::FiniteElementSpace<parallel>> fes_sm;
    //! \brief result_sm is used to transfer data from the result gridfunction
    //! on submesh
    std::shared_ptr<mfem_mgis::GridFunction<parallel>> result_sm;
  };  // end of struct ParaviewExportResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX */
