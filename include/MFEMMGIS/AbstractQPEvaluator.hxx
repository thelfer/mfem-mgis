/*!
 * \file   MFEMMGIS/AbstractQPEvaluator.hxx
 * \brief
 * \author Thomas HElfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEMMGIS_ABSTRACTQPEVALUATOR_HXX
#define LIB_MFEMMGIS_ABSTRACTQPEVALUATOR_HXX

#include <memory>
#include <vector>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;

  /*!
   * \brief structure describing the result of an partial quadrature function
   * evaluator
   *
   * By essence, a partial quadrature function evaluator may want to return a
   * view to an existing partial quadrature function or create a new one.
   */
  struct QPEvaluatorResult : ImmutablePartialQuadratureFunctionView {
    //! \brief constructor to an existing partial quadrature function
    QPEvaluatorResult(const ImmutablePartialQuadratureFunctionView&) noexcept;
    //! \brief constructor from a r-value to a partial quadrature function
    QPEvaluatorResult(PartialQuadratureFunction&&) noexcept;
    //! \brief constructor taking a copy of a partial quadrature function
    explicit QPEvaluatorResult(const PartialQuadratureFunction&) noexcept;
    //! \brief move constructor
    QPEvaluatorResult(QPEvaluatorResult&&) noexcept;
    //! \brief copy constructor
    QPEvaluatorResult(const QPEvaluatorResult&) noexcept;
    //! \brief move assignement
    QPEvaluatorResult& operator=(QPEvaluatorResult&&) noexcept;
    //! \brief destructor
    ~QPEvaluatorResult() noexcept;

   private:
    QPEvaluatorResult& operator=(const QPEvaluatorResult&&) = delete;
    //! \brief internal pointer for memory management, if required
    std::optional<PartialQuadratureFunction> f;
  };

  /*!
   * \brief generic interface for partial quadrature function evaluators
   */
  struct MFEM_MGIS_EXPORT AbstractQPEvaluator {
    [[nodiscard]] virtual const PartialQuadratureSpace& getQuadratureSpace()
        const noexcept = 0;
    [[nodiscard]] virtual std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept = 0;
    //! \return the number of components
    [[nodiscard]] virtual size_type getNumberOfComponents() const noexcept = 0;
    [[nodiscard]] virtual bool isUniform() const noexcept = 0;
    /*!
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx, const real, const real) const noexcept = 0;
    /*!
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<QPEvaluatorResult> evaluate(
        Context&, const real, const real) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractQPEvaluator() noexcept;
  };  // end of AbstractQPEvaluator

  //! \brief options passed to the evaluate function
  struct QPEvaluationOptions {
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
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<QPEvaluatorResult> evaluate(
      Context&,
      const AbstractQPEvaluator&,
      const real,
      const real,
      const QPEvaluationOptions& = QPEvaluationOptions{}) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTQPEVALUATOR_HXX */
