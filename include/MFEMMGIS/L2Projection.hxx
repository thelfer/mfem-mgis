/*!
 * \file   MFEMMGIS/L2Projection.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/01/2026
 */

#ifndef LIB_MFEMMGIS_L2PROJECTION_HXX
#define LIB_MFEMMGIS_L2PROJECTION_HXX

#include <vector>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declarations
  struct LinearSolverHandler;
  struct ImmutablePartialQuadratureFunctionView;

  template <bool parallel>
  struct L2ProjectionResult {};

  /*!
   * \brief compute the L2 projection of the given partial quadrature fields on
   * nodes.
   *
   * \param[in] fcts: functions to be projected
   */
  template <bool parallel>
  std::optional<L2ProjectionResult<parallel>> computeL2Projection(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

  //
  template <>
  MFEM_MGIS_EXPORT std::optional<L2ProjectionResult<true>>
  computeL2Projection<true>(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);
  template <>
  MFEM_MGIS_EXPORT std::optional<L2ProjectionResult<false>>
  computeL2Projection<false>(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_L2PROJECTION_HXX */
