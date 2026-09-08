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

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct MeshDiscretization;
  struct PhysicalSystem;

  //! \brief class meant to extract the values of a grid function along a curve.
  struct MFEM_MGIS_EXPORT PointsSetCurves {
    //! \return a description of each parameters of this class
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] ps: physical system
     * \param[in] params: parameters
     */
    PointsSetCurves(const PhysicalSystem &, const Parameters &);
    /*!
     * \brief constructor
     * \param[in] m: mesh discretization
     * \param[in] params: parameters
     */
    PointsSetCurves(const MeshDiscretization &, const Parameters &);
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
    //! \brief return if the space dimension
    [[nodiscard]] size_type getSpaceDimension() const noexcept;
    //! \brief return if the curvilinear abscissa is calculated
    [[nodiscard]] bool exportCurvilinearAbscissa() const noexcept;
    //! \brief return the curvilinear abscissa
    [[nodiscard]] OptionalReference<const std::vector<real>>
    getCurvilinearAbscissa(Context &) const noexcept;
    //! \brief return if the line curve exports the coordinates
    [[nodiscard]] bool exportCoordinates() const noexcept;
    //! \brief return the coordinates
    [[nodiscard]] std::optional<std::vector<std::vector<real>>> getCoordinates(
        Context &) const noexcept;
    //! \brief return the description of the selected evaluators
    [[nodiscard]] std::vector<std::string> getValuesDescription()
        const noexcept;
    //! \brief return the nodal values of the selected evaluators
    [[nodiscard]] std::optional<std::vector<std::vector<real>>> getValues(
        Context &, const TimeStepStage) const noexcept;
    // \brief destructor
    ~PointsSetCurves() noexcept;

   private:
#ifdef MFEM_USE_MPI
    std::vector<std::pair<
        std::string,
        std::variant<const GridFunction<true> *, const GridFunction<false> *>>>
        gridfunctions;
#else  /* MFEM_USE_MPI */
    std::vector<std::pair<std::string, const GridFunction<false> *>>
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
     * calculated
     */
    bool shallExportCurvilinearAbscissa = true;
    //! \brief flag stating if the coordinates along the line can be retrieved
    bool shallExportCoordinates = false;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVES_HXX */
