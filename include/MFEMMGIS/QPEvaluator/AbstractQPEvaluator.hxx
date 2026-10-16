/*!
 * \file   MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx
 * \brief  This file declares the `AbstractQPEvaluator` class
 * \author Thomas Helfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATOR_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATOR_HXX

#include <memory>
#include <vector>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;

  /*!
   * \brief structure storing a partial quadrature function, if needed
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorResultStorage {
   protected:
    //! \brief default constructor
    QPEvaluatorResultStorage() noexcept;
    /*!
     * \brief move constructor
     * \param[in] src: source
     */
    QPEvaluatorResultStorage(QPEvaluatorResultStorage&& src) noexcept;
    /*!
     * \brief constructor from a partial quadrature function
     * \param[in] f: partial quadrature function
     */
    QPEvaluatorResultStorage(PartialQuadratureFunction&& f) noexcept;
    //
    QPEvaluatorResultStorage(const QPEvaluatorResultStorage&) = delete;
    QPEvaluatorResultStorage& operator=(QPEvaluatorResultStorage&&) = delete;
    QPEvaluatorResultStorage& operator=(const QPEvaluatorResultStorage&) =
        delete;
    //! \brief internal storage for memory management, if required
    std::optional<PartialQuadratureFunction> f;
  };

  /*!
   * \brief structure describing the result of a partial quadrature function
   * evaluator
   *
   * By essence, a partial quadrature function evaluator may want to return a
   * view to an existing partial quadrature function or create a new one.
   *
   * The functions `fromFunction` and `fromView` were added after
   * Issue 360 that showed that the direct usage of the class constructors can
   * be easily misused due to the implicit conversion of a function into an
   * immutable view, see https://github.com/thelfer/mfem-mgis/pull/361 for
   * details
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorResult
      : QPEvaluatorResultStorage,
        ImmutablePartialQuadratureFunctionView {
    /*!
     * \brief build from a function
     * \param[in, out] v: function, moved into the result
     * \return the result
     */
    [[nodiscard]] static QPEvaluatorResult fromFunction(
        PartialQuadratureFunction&& v) noexcept;
    /*!
     * \brief build from a view
     * \param[in] v: view
     * \return the result
     */
    [[nodiscard]] static QPEvaluatorResult fromView(
        const ImmutablePartialQuadratureFunctionView& v) noexcept;
    //! \brief move constructor
    QPEvaluatorResult(QPEvaluatorResult&&) noexcept;
    //
    QPEvaluatorResult& operator=(QPEvaluatorResult&&) noexcept = delete;
    QPEvaluatorResult& operator=(const QPEvaluatorResult&) = delete;
    //! \brief destructor
    ~QPEvaluatorResult() noexcept;

   private:
    //! \brief constructor to an existing partial quadrature function
    explicit QPEvaluatorResult(
        const ImmutablePartialQuadratureFunctionView&) noexcept;
    //! \brief constructor from a r-value to a partial quadrature function
    explicit QPEvaluatorResult(PartialQuadratureFunction&&) noexcept;
  };

  /*!
   * \brief generic interface for partial quadrature function evaluators
   */
  struct MFEM_MGIS_EXPORT AbstractQPEvaluator {
    //! \return the partial quadrature space
    [[nodiscard]] virtual const PartialQuadratureSpace&
    getPartialQuadratureSpace() const noexcept = 0;
    //! \return the number of components
    [[nodiscard]] virtual size_type getNumberOfComponents() const noexcept = 0;
    //! \return if the evaluator is spatially uniform
    [[nodiscard]] virtual bool isUniform() const noexcept = 0;
    /*!
     * \brief return the value of a uniform evaluator
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx,
                    const real t,
                    const real dt) const noexcept = 0;
    /*!
     * \brief evaluate the values at the integration points
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<QPEvaluatorResult> evaluate(
        Context& ctx, const real t, const real dt) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractQPEvaluator() noexcept;
  };  // end of AbstractQPEvaluator

  //! \brief options passed to the evaluate function
  struct QPEvaluationOptions {
    //! \brief if the evaluated values must be stored contiguously
    const bool ensure_contiguous_storage = false;
  };  // end of struct QPEvaluationOptions

  /*!
   * \brief evaluate a partial quadrature function
   *
   * \param[in, out] ctx: execution context
   * \param[in] e: evaluator
   * \param[in] t: time at the beginning of the time step
   * \param[in] dt: time increment
   * \param[in] opts: options
   * \return the result of the evaluation
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]]  //
  std::optional<QPEvaluatorResult>
  evaluate(Context& ctx,
           const AbstractQPEvaluator& e,
           const real t,
           const real dt,
           const QPEvaluationOptions& opts = QPEvaluationOptions{}) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_ABSTRACTQPEVALUATOR_HXX */
