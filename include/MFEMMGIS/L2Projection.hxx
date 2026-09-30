/*!
 * \file   MFEMMGIS/L2Projection.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/01/2026
 */

#ifndef LIB_MFEMMGIS_L2PROJECTION_HXX
#define LIB_MFEMMGIS_L2PROJECTION_HXX

#include <vector>
#include <memory>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declarations
  struct LinearSolverHandler;
  struct ImmutablePartialQuadratureFunctionView;

  //! \brief result of an L2 projection
  template <bool parallel>
  struct L2ProjectionResult {
    /*!
     * \brief submesh created for the resolution. May be empty if the projection
     * is done on the whole mesh.
     */
    std::shared_ptr<SubMesh<parallel>> submesh;
    //! \brief grid function resulting from the projection
    std::unique_ptr<GridFunction<parallel>> result;
  };

  /*!
   * \brief allocate the L2 projection of the given functions
   * \param[in, out] ctx: execution context
   * \param[in] fcts: functions to be projected
   * \return the allocated result
   */
  template <bool parallel>
  [[nodiscard]] std::optional<L2ProjectionResult<parallel>>
  createL2ProjectionResult(
      Context& ctx,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) noexcept;
  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<true>>
  createL2ProjectionResult(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //! \brief sequential specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<false>>
  createL2ProjectionResult(
      Context&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  /*!
   * \brief update the L2 projection of the given partial quadrature fields
   * on nodes.
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] r: result of the projection
   * \param[in, out] s: linear solver handler
   * \param[in] fcts: functions to be projected
   * \return true on success
   */
  template <bool parallel>
  [[nodiscard]] bool updateL2Projection(
      Context& ctx,
      L2ProjectionResult<parallel>& r,
      LinearSolverHandler& s,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) noexcept;
  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool updateL2Projection<true>(
      Context&,
      L2ProjectionResult<true>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //! \brief sequential specialisation
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
   * \param[in, out] ctx: execution context
   * \param[in, out] s: linear solver handler
   * \param[in] fcts: functions to be projected
   * \return the result of the projection
   */
  template <bool parallel>
  std::optional<L2ProjectionResult<parallel>> computeL2Projection(
      Context& ctx,
      LinearSolverHandler& s,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts) noexcept;
  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<L2ProjectionResult<true>>
  computeL2Projection<true>(
      Context&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&) noexcept;
  //! \brief sequential specialisation
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
   * See Peerlings et al. for details.
   *
   * Peerlings, R.H.J., de Borst, R., Brekelmans, W.A.M. and de Vree, J.H.P.
   * (1996). Gradient-Enhanced Damage for Quasi-brittle Materials, International
   * Journal for Numerical Methods in Engineering, 39: 3391-3403.
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] r: result of the regularization
   * \param[in, out] s: linear solver handler
   * \param[in] fcts: functions to be regularized
   * \param[in] l: characteristic length
   * \return true on success
   */
  template <bool parallel>
  [[nodiscard]] bool updateImplicitGradientRegularization(
      Context& ctx,
      ImplicitGradientRegularizationResult<parallel>& r,
      LinearSolverHandler& s,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real l) noexcept;
  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] bool
  updateImplicitGradientRegularization<true>(
      Context&,
      ImplicitGradientRegularizationResult<true>&,
      LinearSolverHandler&,
      const std::vector<ImmutablePartialQuadratureFunctionView>&,
      const real) noexcept;
  //! \brief sequential specialisation
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
   * \param[in, out] ctx: execution context
   * \param[in, out] s: linear solver handler
   * \param[in] fcts: functions to be regularized
   * \param[in] l: characteristic length
   * \return the result of the regularization
   */
  template <bool parallel>
  std::optional<ImplicitGradientRegularizationResult<parallel>>
  computeImplicitGradientRegularization(
      Context& ctx,
      LinearSolverHandler& s,
      const std::vector<ImmutablePartialQuadratureFunctionView>& fcts,
      const real l) noexcept;
  //! \brief parallel specialisation
  template <>
  MFEM_MGIS_EXPORT
      [[nodiscard]] std::optional<ImplicitGradientRegularizationResult<true>>
      computeImplicitGradientRegularization<true>(
          Context&,
          LinearSolverHandler&,
          const std::vector<ImmutablePartialQuadratureFunctionView>&,
          const real) noexcept;
  //! \brief sequential specialisation
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
