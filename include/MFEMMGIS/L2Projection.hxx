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
  struct L2ProjectionResult {
    /*!
     * \brief submesh created for the resolution. May be empty if the projection
     * is done on the whole mesh.
     */
    std::unique_ptr<SubMesh<parallel>> submesh;
    //! \brief finite element space.
    std::unique_ptr<FiniteElementSpace<parallel>> fe_space;
    //! \brief grid function resulting from the projection
    std::unique_ptr<GridFunction<parallel>> result;
  };

  template <bool parallel>
  [[nodiscard]] std::optional<L2ProjectionResult<parallel>>
  createL2ProjectionResult(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<true>>
  createL2ProjectionResult(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<false>>
  createL2ProjectionResult(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  /*!
   * \brief update the L2 projection of the given partial quadrature fields
   * on nodes.
   *
   * \param[in] fcts: functions to be projected
   */
  template <bool parallel>
  [[nodiscard]] bool updateL2Projection(
      Context&,
      L2ProjectionResult<parallel>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool updateL2Projection<true>(
      Context&,
      L2ProjectionResult<true>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool updateL2Projection<false>(
      Context&,
      L2ProjectionResult<false>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  /*!
   * \brief compute the L2 projection of the given partial quadrature fields
   * on nodes.
   * This function first calls `createL2ProjectionResult` and then
   * `updateL2Projection`.
   *
   * \param[in] fcts: functions to be projected
   */
  template <bool parallel>
  std::optional<L2ProjectionResult<parallel>> computeL2Projection(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<true>>
  computeL2Projection<true>(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<false>>
  computeL2Projection<false>(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_L2PROJECTION_HXX */
