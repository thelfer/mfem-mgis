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

  //! \brief a simple alias
  template <size_type N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      using Point = ::tfel::math::tvector<N, real>;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context& ctx, const Parameter& p) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context& ctx, const Parameter& p) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context& ctx, const Parameters& p) noexcept;

  template <unsigned short N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::string toString(const Point<N>& pt) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context& ctx,
          const std::map<std::string, std::vector<Point<N>>, std::less<>>&
              pointsSets,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameter& p) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context& ctx,
          const std::map<std::string, Point<N>, std::less<>>& pts,
          const Parameters& p) noexcept;

  /*!
   * \return the curvilinear abscissae along the given points set
   * \param[in] pts: points set
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<real> computeCurvilinearAbscissae(
      const std::vector<Point<2>>& pts) noexcept;
  /*!
   * \return the curvilinear abscissae along the given points set
   * \param[in] pts: points set
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<real> computeCurvilinearAbscissae(
      const std::vector<Point<3>>& pts) noexcept;

  // partial specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&, const Parameters&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&, const Parameters&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::string toString<2>(
      const Point<2>&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::string toString<3>(
      const Point<3>&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&,
      const std::map<std::string, Point<2>, std::less<>>&,
      const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(
      Context&,
      const std::map<std::string, std::vector<Point<2>>, std::less<>>&,
      const std::map<std::string, Point<2>, std::less<>>&,
      const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&,
                       const std::map<std::string, Point<2>, std::less<>>&,
                       const Parameters&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&,
      const std::map<std::string, Point<3>, std::less<>>&,
      const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(
      Context&,
      const std::map<std::string, std::vector<Point<3>>, std::less<>>&,
      const std::map<std::string, Point<3>, std::less<>>&,
      const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&,
                       const std::map<std::string, Point<3>, std::less<>>&,
                       const Parameters&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_GEOMETRY_HXX */
