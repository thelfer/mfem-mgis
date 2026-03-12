/*!
 * \file   MFEMMGIS/AbstractPartialQuadratureFunctionEvaluator.hxx
 * \brief
 * \author Thomas HElfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX
#define LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX

#include <memory>
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
  struct PartialQuadratureFunctionEvaluatorResult
      : ImmutablePartialQuadratureFunctionView {
    //! \brief constructor to an existing partial quadrature function
    PartialQuadratureFunctionEvaluatorResult(
        const ImmutablePartialQuadratureFunctionView&) noexcept;
    //! \brief constructor from a r-value to a partial quadrature function
    PartialQuadratureFunctionEvaluatorResult(
        PartialQuadratureFunction&&) noexcept;
    //! \brief constructor taking a copy of a partial quadrature function
    explicit PartialQuadratureFunctionEvaluatorResult(
        const PartialQuadratureFunction&) noexcept;
    //! \brief move constructor
    PartialQuadratureFunctionEvaluatorResult(
        PartialQuadratureFunctionEvaluatorResult&&) noexcept;
    //! \brief copy constructor
    PartialQuadratureFunctionEvaluatorResult(
        const PartialQuadratureFunctionEvaluatorResult&) noexcept;
    //! \brief move assignement
    PartialQuadratureFunctionEvaluatorResult& operator=(
        PartialQuadratureFunctionEvaluatorResult&&) noexcept;
    //! \brief destructor
    ~PartialQuadratureFunctionEvaluatorResult() noexcept;

   private:
    PartialQuadratureFunctionEvaluatorResult& operator=(
        const PartialQuadratureFunctionEvaluatorResult&&) = delete;
    //! \brief internal pointer for memory management, if required
    std::optional<PartialQuadratureFunction> f;
  };

  /*!
   * \brief generic interface for partial quadrature function evaluators
   */
  struct MFEM_MGIS_EXPORT AbstractPartialQuadratureFunctionEvaluator {
    [[nodiscard]] virtual const PartialQuadratureSpace& getQuadratureSpace()
        const noexcept = 0;
    [[nodiscard]] virtual std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept = 0;
    //! \return the number of components
    [[nodiscard]] virtual size_type getNumberOfComponents() const noexcept = 0;
    //!
    [[nodiscard]] virtual std::optional<
        PartialQuadratureFunctionEvaluatorResult>
    evaluate(Context&) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractPartialQuadratureFunctionEvaluator() noexcept;
  };  // end of AbstractPartialQuadratureFunctionEvaluator

  //! \brief options passed to the evaluate function
  struct PartialQuadratureFunctionEvaluationOptions {
    const bool ensure_contiguous_storage = false;
  };  // end of struct PartialQuadratureFunctionEvaluationOptions

  /*!
   * \brief evaluate a partial quadrature function
   *
   * \param[in, out] ctx: execution context
   * \param[in] e: evaluator
   * \param[in] opts: options
   */
  MFEM_MGIS_EXPORT
  [[nodiscard]]  //
  std::optional<PartialQuadratureFunctionEvaluatorResult>
  evaluate(Context&,
           const AbstractPartialQuadratureFunctionEvaluator&,
           const PartialQuadratureFunctionEvaluationOptions& =
               PartialQuadratureFunctionEvaluationOptions{}) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX */
