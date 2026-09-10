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
   * \brief structure in charge of interpolateing values at
   * given points, mostly for post-processing purposes.
   *
   * This structure uses features provided by the gslib libray and MFEM shall be
   * built around it.
   */
  struct MFEM_MGIS_EXPORT GridFunctionInterpolator {
    /*!
     * \brief default constructor from a finite element discretization
     * \param[in] m: finite element spaces manager
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    explicit GridFunctionInterpolator(
        const FiniteElementSpacesManager&) noexcept;
    /*!
     * \brief constructor from a set of 2D points
     *
     * \param[in] ctx: execution context
     * \param[in] m: finite element spaces manager
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context&,
                             const FiniteElementSpacesManager&,
                             const std::vector<Point<2>>&);
    /*!
     * \brief constructor from a set of 3D points
     *
     * \param[in] ctx: execution context
     * \param[in] m: finite element spaces manager
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context&,
                             const FiniteElementSpacesManager&,
                             const std::vector<Point<3>>&);
    /*!
     * \brief default constructor from a finite element discretization
     * \param[in] fed: finite element discretization
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    explicit GridFunctionInterpolator(
        const FiniteElementDiscretization&) noexcept;
    /*!
     * \brief constructor from a set of 2D points
     *
     * \param[in] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context&,
                             const FiniteElementDiscretization&,
                             const std::vector<Point<2>>&);
    /*!
     * \brief constructor from a set of 3D points
     *
     * \param[in] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] pts: points to be added
     *
     * \note this constructor is provided for simplying the declaration
     * of an interpolator as only the underlying finite element spaces manager
     * is required
     */
    GridFunctionInterpolator(Context&,
                             const FiniteElementDiscretization&,
                             const std::vector<Point<3>>&);
    /*!
     * \brief add the given 2D points to the list of points to be post-processed
     *
     * \param[in] ctx: execution context
     * \param[in] pts: points to be added
     */
    [[nodiscard]] bool addPoints(Context&,
                                 const std::vector<Point<2>>&) noexcept;
    /*!
     * \brief add the given 3D points to the list of points to be post-processed
     *
     * \param[in] ctx: execution context
     * \param[in] pts: points to be added
     */
    [[nodiscard]] bool addPoints(Context&,
                                 const std::vector<Point<3>>&) noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \brief interpolate the given grid function at the previously defined
     * points
     *
     * \param[in] ctx: execution context
     * \param[in] f: function to be interpolated
     *
     * \note points are searched each time since the we have to call
     * `SetNodalFESpace` on the underlying mesh and that this mesh may also be
     * used to build order finite element spaces.
     */
    [[nodiscard]] std::optional<tfel::math::matrix<real>> interpolate(
        Context&, const GridFunction<true>&) noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief interpolate the given grid function at the previously defined
     * points
     *
     * \param[in] ctx: execution context
     * \param[in] f: function to be interpolated
     *
     * \note points are searched each time since the we have to call
     * `SetNodalFESpace` on the underlying mesh and that this mesh may also be
     * used to build order finite element spaces.
     */
    [[nodiscard]] std::optional<tfel::math::matrix<real>> interpolate(
        Context&, const GridFunction<false>&) noexcept;

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
