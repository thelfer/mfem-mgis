/*!
 * \file   MFEMMGIS/PostProcessing/PointsSetCurvesPostProcessing.hxx
 * \brief  This file declares the `PointsSetCurvesPostProcessing` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESPOSTPROCESSING_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESPOSTPROCESSING_HXX

#include <fstream>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx"
#include "MFEMMGIS/PostProcessing/PostProcessingBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;
  struct Context;

  /*!
   * \brief post-processing meant to export the values of grid functions on a
   * points set to a file.
   */
  struct MFEM_MGIS_EXPORT PointsSetCurvesPostProcessing
      : public PostProcessingBase {
    //! \return a description of each parameter of this struct
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \return a description of the post-processing
    static std::string getDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] ps: physical system
     * \param[in] manager: finite element spaces manager
     * \param[in] params: parameters
     */
    PointsSetCurvesPostProcessing(PhysicalSystem& ps,
                                  const FiniteElementSpacesManager& manager,
                                  const Parameters& params);
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
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context& ctx, const real t) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(
        Context& ctx,
        const TimeStep& ts,
        const bool isPostProcessingRequired) noexcept override;
    //! \brief destructor
    ~PointsSetCurvesPostProcessing() noexcept override;

   private:
    //! \brief curve writer
    PointsSetCurvesWriter writer;
    /*!
     * \brief boolean stating if values shall be exported at the beginning of
     * the first time step
     */
    const bool executeInitialPostProcessing;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESPOSTPROCESSING_HXX */
