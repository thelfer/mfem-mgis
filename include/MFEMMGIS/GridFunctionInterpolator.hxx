/*!
 * \file   MFEMMGIS/GridFunctionInterpolator.hxx
 * \brief  This file declares the `GridFunctionInterpolator` class
 * \author Thomas Helfer
 * \date   04/09/2026
 */

#ifndef LIB_MFEMMGIS_GRIDFUNCTIONINTERPOLATOR_HXX
#define LIB_MFEMMGIS_GRIDFUNCTIONINTERPOLATOR_HXX

#ifndef MFEM_USE_GSLIB
#error "gslib support in MFEM is not enabled"
#endif

#ifndef MGIS_HAVE_TFEL
#error "TFEL support in MGIS is not enabled"
#endif

#include <vector>
#include <optional>
#include "TFEL/Math/matrix.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Geometry.hxx"
#include "MFEMMGIS/FiniteElementSpacesManager.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief structure in charge of interpolating values at
   * given points, mostly for post-processing purposes.
   *
   * This structure uses features provided by the gslib library and MFEM shall
   * be built with it.
   */
  struct MFEM_MGIS_EXPORT GridFunctionInterpolator {
    /*!
     * \brief constructor from a finite element spaces manager
     * \param[in] m: finite element spaces manager
     */
    explicit GridFunctionInterpolator(
        const FiniteElementSpacesManager& m) noexcept;
    /*!
     * \brief constructor from a set of 2D points
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: finite element spaces manager
     * \param[in] pts: points to be added
     */
    GridFunctionInterpolator(Context& ctx,
                             const FiniteElementSpacesManager& m,
                             const std::vector<Point<2>>& pts);
    /*!
     * \brief constructor from a set of 3D points
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: finite element spaces manager
     * \param[in] pts: points to be added
     */
    GridFunctionInterpolator(Context& ctx,
                             const FiniteElementSpacesManager& m,
                             const std::vector<Point<3>>& pts);
    /*!
     * \brief constructor from a finite element discretization
     * \param[in] fed: finite element discretization
     *
     * \note this constructor is provided for simplifying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    explicit GridFunctionInterpolator(
        const FiniteElementDiscretization& fed) noexcept;
    /*!
     * \brief constructor from a set of 2D points
     *
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplifying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context& ctx,
                             const FiniteElementDiscretization& fed,
                             const std::vector<Point<2>>& pts);
    /*!
     * \brief constructor from a set of 3D points
     *
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplifying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context& ctx,
                             const FiniteElementDiscretization& fed,
                             const std::vector<Point<3>>& pts);
    /*!
     * \brief add the given 2D points to the list of points to be post-processed
     *
     * \param[in, out] ctx: execution context
     * \param[in] pts: points to be added
     * \return true on success
     */
    [[nodiscard]] bool addPoints(Context& ctx,
                                 const std::vector<Point<2>>& pts) noexcept;
    /*!
     * \brief add the given 3D points to the list of points to be post-processed
     *
     * \param[in, out] ctx: execution context
     * \param[in] pts: points to be added
     * \return true on success
     */
    [[nodiscard]] bool addPoints(Context& ctx,
                                 const std::vector<Point<3>>& pts) noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \brief interpolate the given grid function at the previously defined
     * points
     *
     * \param[in, out] ctx: execution context
     * \param[in] f: function to be interpolated
     * \return the interpolated values, one row per point and one column per
     * component
     *
     * \note points are searched each time. The underlying mesh is given
     * nodes of order 1 if it has none (`EnsureNodes`), which keeps its geometry
     * and the integration rules of the behaviour integrators.
     */
    [[nodiscard]] std::optional<tfel::math::matrix<real>> interpolate(
        Context& ctx, const GridFunction<true>& f) noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief interpolate the given grid function at the previously defined
     * points
     *
     * \param[in, out] ctx: execution context
     * \param[in] f: function to be interpolated
     * \return the interpolated values, one row per point and one column per
     * component
     *
     * \note points are searched each time. The underlying mesh is given
     * nodes of order 1 if it has none (`EnsureNodes`), which keeps its geometry
     * and the integration rules of the behaviour integrators.
     */
    [[nodiscard]] std::optional<tfel::math::matrix<real>> interpolate(
        Context& ctx, const GridFunction<false>& f) noexcept;

    //! \brief destructor
    ~GridFunctionInterpolator();

   private:
    //! \brief underlying finite element space manager
    FiniteElementSpacesManager fespaces_manager;
    //! \brief list of points stored byVDIM (XYXY... in 2D, XYZXYZ.. in 3D)
    std::vector<real> points;
    //!
  };  // end of struct GridFunctionInterpolator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_GRIDFUNCTIONINTERPOLATOR_HXX */
