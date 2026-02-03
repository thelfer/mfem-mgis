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

  //! \brief a simple alias
  template <bool parallel>
  using ImplicitGradientRegularizationResult = L2ProjectionResult<parallel>;

  /*!
   * \brief update the implicit gradient regularisation of the given partial
   * quadrature fields on nodes.
   *
   * For a scalar function \f$f\f$, the implicit gradient regularization
   * \f$\bar{f}\f$
   * is defined as the solution of:
   *
   * \f[
   * \bar{f}-l^{2}\cdot\Delta\bar{f}=f
   * \f]
   *
   * where \f$l\f$ is a characteristic length.
   *
   * See Peerling et al. for details.
   *
   * Peerlings, R.H.J., de Borst, R., Brekelmans, W.A.M. and de Vree, J.H.P.
   * (1996). Gradient-Enhanced Damage for Quasi-brittle Materials, International
   * Journal for Numerical Methods in Engineering, 39: 3391 3403.
   *
   * \param[in] ctx: execution context
   * \param[in] l: linear solver handler
   * \param[in] fcts: functions to be reguarlized
   * \param[in] l: characteristic length
   */
  template <bool parallel>
  [[nodiscard]] bool updateImplicitGradientRegularization(
      Context&,
      ImplicitGradientRegularizationResult<parallel>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const real) noexcept;
  //
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  updateImplicitGradientRegularization<true>(
      Context&,
      ImplicitGradientRegularizationResult<true>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const real) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  updateImplicitGradientRegularization<false>(
      Context&,
      ImplicitGradientRegularizationResult<false>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const real) noexcept;
  /*!
   * \brief compute the implicit gradient regularisation of the given partial
   * quadrature fields on nodes. This function first calls
   * `createL2ProjectionResult` and then `updateImplicitGradientRegularization`.
   *
   * \param[in] ctx: execution context
   * \param[in] l: linear solver handler
   * \param[in] fcts: functions to be reguarlized
   * \param[in] l: characteristic length
   */
  template <bool parallel>
  std::optional<ImplicitGradientRegularizationResult<parallel>>
  computeImplicitGradientRegularization(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const real) noexcept;
  //
  template <>
  MFEM_MGIS_EXPORT
      [[nodiscard]] std::optional<ImplicitGradientRegularizationResult<true>>
      computeImplicitGradientRegularization<true>(
          Context&,
          LinearSolverHandler&,
          const std::vector<ImmutablePartialQuadratureFunctionView>&,
          const real) noexcept;
  template <>
  MFEM_MGIS_EXPORT
      [[nodiscard]] std::optional<ImplicitGradientRegularizationResult<false>>
      computeImplicitGradientRegularization<false>(
          Context&,
          LinearSolverHandler&,
          const std::vector<ImmutablePartialQuadratureFunctionView>&,
          const real) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_L2PROJECTION_HXX */
