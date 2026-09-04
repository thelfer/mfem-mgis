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
          Context&, const Parameter&) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context&, const Parameters&) noexcept;

  // partial specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&, const Parameters&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&, const Parameters&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_GEOMETRY_HXX */
