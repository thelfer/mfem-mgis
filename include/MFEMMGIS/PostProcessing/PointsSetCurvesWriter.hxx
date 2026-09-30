/*!
 * \file   MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx
 * \brief  This file declares the `PointsSetCurvesWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX

#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"
#include "MFEMMGIS/Utilities/DataFileUtilities.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurves.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct FiniteElementDiscretization;

  /*!
   * \brief helper class meant to write the values of grid functions on a
   * points set to an output file
   */
  struct MFEM_MGIS_EXPORT PointsSetCurvesWriter {
    //! \return a description of each parameter
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] manager: finite element spaces manager
     * \param[in] parameters: parameters
     */
    PointsSetCurvesWriter(const FiniteElementSpacesManager& manager,
                          const Parameters& parameters);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] parameters: parameters
     */
    PointsSetCurvesWriter(const FiniteElementDiscretization& fed,
                          const Parameters& parameters);
#ifdef MFEM_USE_MPI
    /*!
     * \brief add a grid function (parallel version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           const GridFunction<true>& f) noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add a grid function (sequential version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           const GridFunction<false>& f) noexcept;
    /*!
     * \brief write the file headers
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    bool writeFileHeader(Context& ctx);
    /*!
     * \brief write the values
     * \param[in, out] ctx: execution context
     * \param[in] ts: time step
     * \param[in] tss: time step stage
     * \return true on success
     */
    bool writeValues(Context& ctx,
                     const TimeStep& ts,
                     const TimeStepStage& tss);

   private:
    //! \brief underlying finite element space manager
    FiniteElementSpacesManager fespaces_manager;
    //! \brief points set curves
    PointsSetCurves curves;
    //! \brief output file
    std::ofstream out;
    //! \brief boolean stating if new curves are allowed
    bool allowNewPointsSetCurves = true;
  };  // end of PointsSetCurvesWriter

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX */
