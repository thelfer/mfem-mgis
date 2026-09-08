/*!
 * \file   MFEMMGIS/PostProcessing/PointsSetCurvesPostProcessing.hxx
 * \brief  This file declares the `PointsSetCurvesPostProcessing` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_CURVESPOSTPROCESSING_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_CURVESPOSTPROCESSING_HXX

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
   * \brief post-processing meant to export values extracted using
   * instances of the `AbstractCurve` struct to a file.
   */
  struct MFEM_MGIS_EXPORT PointsSetCurvesPostProcessing
      : public PostProcessingBase {
    //! \return a description of each parameters of this struct
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \return a description of the post-processing
    static std::string getDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] params: parameters
     */
    PointsSetCurvesPostProcessing(PhysicalSystem &, const Parameters &);
#ifdef MFEM_USE_MPI
    /*!
     * \brief add a grid function  (parallel version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     */
    [[nodiscard]] bool add(Context &,
                           std::string_view,
                           const GridFunction<true> &) noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add a grid function (sequential version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     */
    [[nodiscard]] bool add(Context &,
                           std::string_view,
                           const GridFunction<false> &) noexcept;
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &, const real) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(Context &,
                                                  const TimeStep &,
                                                  const bool) noexcept override;
    // \brief destructor
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

#endif /* LIB_MFEMMGIS_POSTPROCESSING_CURVESPOSTPROCESSING_HXX */
