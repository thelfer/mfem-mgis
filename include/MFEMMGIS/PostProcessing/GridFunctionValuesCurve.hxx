/*!
 * \file   MFEMMGIS/PostProcessing/GridFunctionValuesCurves.hxx
 * \brief  This file declares the `GridFunctionValuesCurve` class.
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_GRIDFUNCTIONVALUESCURVES_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_GRIDFUNCTIONVALUESCURVES_HXX

#include <map>
#include <string>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
#include "MFEMMGIS/Geometry.hxx"
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
#include "MFEMMGIS/PostProcessing/AbstractCurve.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct PhysicalSystem;

  /*!
   * \brief curve retrieving values of a grid function at some points.
   *
   * This class is mostly a wrapper around the GridFunctionInterpolator class
   */
  struct MFEM_MGIS_EXPORT GridFunctionValuesCurve : public AbstractCurve {
    //! \return a description of each parameters
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \return a description of the curve
    [[nodiscard]] static std::string getDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] params: parameters
     */
    GridFunctionValuesCurve(Context &, PhysicalSystem &, const Parameters &);
    /*!
     * \brief set the grid function to be interpolated (parallel case)
     * \param[in] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] fct: grid function
     */
    [[nodiscard]] virtual bool setGridRunction(
        Context &, std::string_view, const GridFunction<true> &) noexcept;
    /*!
     * \brief set the grid function to be interpolated (sequential case)
     * \param[in] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] fct: grid function
     */
    [[nodiscard]] virtual bool setGridRunction(
        Context &, std::string_view, const GridFunction<false> &) noexcept;
    /*!
     * \brief add points ont which the grid function is to be interpolated
     * \param[in] ctx: execution context
     * \param[in] parameter: parameter defining the set of points
     */
    [[nodiscard]] virtual bool addPoints(Context &, const Parameter &) noexcept;
    //
    [[nodiscard]] std::vector<std::string> getDescriptions()
        const noexcept override;
    [[nodiscard]] std::optional<std::vector<real>> getValues(
        Context &ctx, const TimeStepStage) const noexcept override;
    //! \brief destructor
    ~GridFunctionValuesCurve() noexcept override;

   private:
    //! \return if the points are defined
    [[nodiscard]] bool arePointsDefined() const noexcept;
    //! \brief underlying physical system
    PhysicalSystem &physicalSystem;
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    //! \brief list of points
    std::variant<std::vector<Point<2>>, std::vector<Point<3>>> points;
#endif /* MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR */
    //! \brief name of the grid function
    std::string name;
    /*!
     * \brief pointer to the grid function from which values are interpolated
     * (parallel case)
     */
    const GridFunction<true> *parallel_fct = nullptr;
    /*!
     * \brief pointer to the grid function from which values are interpolated
     * (sequential case)
     */
    const GridFunction<false> *sequential_fct = nullptr;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_GRIDFUNCTIONVALUESCURVES_HXX */
