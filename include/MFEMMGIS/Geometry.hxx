/*!
 * \file   MFEMMGIS/Geometry.hxx
 * \brief  This header declares some utility functions to declare points and
 * lines from `Parameters`
 *
 * \author Thomas Helfer
 * \date   04/09/2026
 */

#ifndef LIB_MFEMMGIS_GEOMETRY_HXX
#define LIB_MFEMMGIS_GEOMETRY_HXX

#ifndef MGIS_HAVE_TFEL
#error "TFEL support in MGIS is required"
#endif /* MGIS_HAVE_TFEL */

#include <map>
#include <string>
#include <vector>
#include "TFEL/Math/tvector.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  //! \brief a point in dimension N
  template <size_type N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      using Point = ::tfel::math::tvector<N, real>;

  /*!
   * \brief create a point from a parameter
   * \return the point defined by the given parameter, empty on failure
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter holding the coordinates of the point
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context& ctx, const Parameter& p) noexcept;

  /*!
   * \brief create a points set from a parameter
   * \return the points set defined by the given parameter, empty on failure
   * \param[in, out] ctx: execution context
   * \param[in] p: list of points or parameters defining a curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context& ctx, const Parameter& p) noexcept;

  /*!
   * \brief discretize a curve
   * \return the points of the curve defined by the given parameters, empty on
   * failure
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters defining the curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context& ctx, const Parameters& p) noexcept;

  /*!
   * \brief convert a point to a string
   * \return a string representation of the given point
   * \param[in] pt: point
   */
  template <unsigned short N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::string toString(const Point<N>& pt) noexcept;

  /*!
   * \brief create a point from a parameter
   * \return the point defined by the given parameter, empty on failure.
   * The parameter may be the name of one of the given points.
   * \param[in, out] ctx: execution context
   * \param[in] pts: named points
   * \param[in] p: name of a point or coordinates of the point
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept;

  /*!
   * \brief create a points set from a parameter
   * \return the points set defined by the given parameter, empty on failure.
   * The parameter may be the name of one of the given points sets.
   * The points may be given by their names.
   * \param[in, out] ctx: execution context
   * \param[in] pointsSets: named points sets
   * \param[in] pts: named points
   * \param[in] p: name of a points set, list of points or parameters defining
   * a curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context& ctx,
          const std::map<std::string, std::vector<Point<N>>, std::less<>>&
              pointsSets,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept;

  /*!
   * \brief discretize a curve
   * \return the points of the curve defined by the given parameters, empty on
   * failure.
   * The points defining the curve may be given by their names.
   * \param[in, out] ctx: execution context
   * \param[in] pts: named points
   * \param[in] p: parameters defining the curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameters& p) noexcept;

  /*!
   * \brief compute the curvilinear abscissae along a points set
   * \return the curvilinear abscissae along the given points set
   * \param[in] pts: points set
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<real> computeCurvilinearAbscissae(
      const std::vector<Point<2>>& pts) noexcept;
  /*!
   * \brief compute the curvilinear abscissae along a points set
   * \return the curvilinear abscissae along the given points set
   * \param[in] pts: points set
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<real> computeCurvilinearAbscissae(
      const std::vector<Point<3>>& pts) noexcept;

  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&, const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&, const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(Context&, const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(Context&, const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&, const Parameters&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&, const Parameters&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::string toString<2>(
      const Point<2>&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::string toString<3>(
      const Point<3>&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&,
      const std::map<std::string, Point<2>, std::less<>>&,
      const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(
      Context&,
      const std::map<std::string, std::vector<Point<2>>, std::less<>>&,
      const std::map<std::string, Point<2>, std::less<>>&,
      const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&,
                       const std::map<std::string, Point<2>, std::less<>>&,
                       const Parameters&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&,
      const std::map<std::string, Point<3>, std::less<>>&,
      const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(
      Context&,
      const std::map<std::string, std::vector<Point<3>>, std::less<>>&,
      const std::map<std::string, Point<3>, std::less<>>&,
      const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&,
                       const std::map<std::string, Point<3>, std::less<>>&,
                       const Parameters&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_GEOMETRY_HXX */
