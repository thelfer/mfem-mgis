/*!
 * \file   MFEMMGIS/PostProcessing/PointsSetCurves.hxx
 * \brief  This file declares the `PointsSetCurves` class
 * \author Thomas Helfer
 * \date   24/03/2025
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVES_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVES_HXX

#include <map>
#include <vector>
#include <variant>
#include <utility>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#ifdef MGIS_HAVE_TFEL
#include "MFEMMGIS/Geometry.hxx"
#endif /* MGIS_HAVE_TFEL */
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct MeshDiscretization;
  struct PhysicalSystem;

  //! \brief class meant to extract the values of grid functions at points.
  struct MFEM_MGIS_EXPORT PointsSetCurves {
    //! \return a description of each parameter of this class
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] manager: finite element spaces manager
     * \param[in] params: parameters
     */
    PointsSetCurves(const FiniteElementSpacesManager& manager,
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
    //! \return the space dimension
    [[nodiscard]] size_type getSpaceDimension() const noexcept;
    //! \return if the curvilinear abscissa shall be exported
    [[nodiscard]] bool exportCurvilinearAbscissa() const noexcept;
    /*!
     * \brief return the curvilinear abscissa
     * \param[in, out] ctx: execution context
     * \return the curvilinear abscissa, empty if not exported
     */
    [[nodiscard]] OptionalReference<const std::vector<real>>
    getCurvilinearAbscissa(Context& ctx) const noexcept;
    //! \return if the line curve exports the coordinates
    [[nodiscard]] bool exportCoordinates() const noexcept;
    /*!
     * \brief return the coordinates of the points
     * \param[in, out] ctx: execution context
     * \return the coordinates of the points, one vector per component
     */
    [[nodiscard]] std::optional<std::vector<std::vector<real>>> getCoordinates(
        Context& ctx) const noexcept;
    //! \return the description of the values of the grid functions
    [[nodiscard]] std::vector<std::string> getValuesDescription()
        const noexcept;
    /*!
     * \brief return the values of the grid functions at the points
     * \param[in, out] ctx: execution context
     * \param[in] ts: time step stage, unused
     * \return the values of the grid functions, one vector per component
     */
    [[nodiscard]] std::optional<std::vector<std::vector<real>>> getValues(
        Context& ctx, const TimeStepStage ts) const noexcept;
    //! \brief destructor
    ~PointsSetCurves() noexcept;

   private:
    //! \brief underlying finite element space manager
    FiniteElementSpacesManager fespaces_manager;
    //! \brief list of registered grid functions
#ifdef MFEM_USE_MPI
    std::vector<std::pair<
        std::string,
        std::variant<const GridFunction<true>*, const GridFunction<false>*>>>
        gridfunctions;
#else  /* MFEM_USE_MPI */
    std::vector<std::pair<std::string, const GridFunction<false>*>>
        gridfunctions;
#endif /* MFEM_USE_MPI */
#ifdef MGIS_HAVE_TFEL
    //! \brief points set
    std::variant<std::vector<Point<2>>, std::vector<Point<3>>> points;
#endif /* MGIS_HAVE_TFEL */
    //! \brief curvilinear abscissae along the line
    std::vector<real> curvilinearAbscissae;
    /*!
     * \brief flag stating if the curvilinear abscissa of the line are
     * exported
     */
    bool shallExportCurvilinearAbscissa = true;
    //! \brief flag stating if the coordinates along the line can be retrieved
    bool shallExportCoordinates = false;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVES_HXX */
