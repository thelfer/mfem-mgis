/*!
 * \file   MFEMMGIS/PostProcessing/GridFunctionValuesCurve.hxx
 * \brief  This file declares the `GridFunctionValuesCurve` class.
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_GRIDFUNCTIONVALUESCURVES_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_GRIDFUNCTIONVALUESCURVES_HXX

#include <map>
#include <string>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
#include "MFEMMGIS/Geometry.hxx"
#include "MFEMMGIS/GridFunctionInterpolator.hxx"
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
    //! \return a description of each parameter
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \return a description of the curve
    [[nodiscard]] static std::string getDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] manager: finite element spaces manager
     * \param[in] parameters: parameters
     */
    GridFunctionValuesCurve(Context &ctx,
                            PhysicalSystem &ps,
                            const FiniteElementSpacesManager &manager,
                            const Parameters &parameters);
    /*!
     * \brief set the grid function to be interpolated (parallel case)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     * \return true on success
     */
    [[nodiscard]] virtual bool setGridRunction(
        Context &ctx, std::string_view n, const GridFunction<true> &f) noexcept;
    /*!
     * \brief set the grid function to be interpolated (sequential case)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     * \return true on success
     */
    [[nodiscard]] virtual bool setGridRunction(
        Context &ctx,
        std::string_view n,
        const GridFunction<false> &f) noexcept;
    /*!
     * \brief set the points on which the grid function is to be interpolated
     * \param[in, out] ctx: execution context
     * \param[in] parameter: parameter defining the set of points
     * \return true on success
     *
     * \note previously defined points are replaced
     */
    [[nodiscard]] virtual bool addPoints(Context &ctx,
                                         const Parameter &parameter) noexcept;
    //
    [[nodiscard]] std::vector<std::string> getDescriptions()
        const noexcept override;
    /*!
     * \brief interpolate the grid function at the points
     * \param[in, out] ctx: execution context
     * \param[in] ts: time step stage, ignored
     * \return the values of the grid function at the points
     */
    [[nodiscard]] std::optional<std::vector<real>> getValues(
        Context &ctx, const TimeStepStage ts) const noexcept override;
    //! \brief destructor
    ~GridFunctionValuesCurve() noexcept override;

   private:
    //! \return if the points are defined
    [[nodiscard]] bool arePointsDefined() const noexcept;
    //! \brief underlying physical system
    PhysicalSystem &physicalSystem;
    //! \brief underlying finite element space manager
    FiniteElementSpacesManager fespaces_manager;
#ifdef MFEMMGIS_HAVE_GSLIBGRIDFUNCTIONINTERPOLATOR
    //! \brief list of points
    std::variant<std::vector<Point<2>>, std::vector<Point<3>>> points;
    //! \brief interpolator, built at the first call to `getValues`
    mutable std::optional<GridFunctionInterpolator> interpolator;
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
