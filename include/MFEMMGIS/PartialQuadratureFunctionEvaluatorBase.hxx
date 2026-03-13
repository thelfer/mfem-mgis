/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionEvaluatorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   12/03/2026
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATORBASE_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATORBASE_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractPartialQuadratureFunctionEvaluator.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for most partial quadrature evaluators
   *
   * \note by default, evaluators are assumed to be spatially variable, thus the
   * default implementation of `isUniform` return false and `getUniformValue`
   * returns an error.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureFunctionEvaluatorBase
      : AbstractPartialQuadratureFunctionEvaluator {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     */
    PartialQuadratureFunctionEvaluatorBase(
        std::shared_ptr<const PartialQuadratureSpace>);
    //
    [[nodiscard]] const PartialQuadratureSpace& getQuadratureSpace()
        const noexcept override;
    [[nodiscard]] std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept override;
    [[nodiscard]] bool isUniform() const noexcept override;
    [[nodiscard]] std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx,
                    const real,
                    const real) const noexcept override;
    //! \brief destructor
    ~PartialQuadratureFunctionEvaluatorBase() noexcept override;

   protected:
    //! \brief partial quadrature space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
  };  // end of struct PartialQuadratureFunctionEvaluatorBase

  //! \brief base class for evaluators returning uniform scalar values
  struct MFEM_MGIS_EXPORT UniformScalarPartialQuadratureFunctionEvaluatorBase
      : PartialQuadratureFunctionEvaluatorBase {
    using PartialQuadratureFunctionEvaluatorBase::
        PartialQuadratureFunctionEvaluatorBase;
    //
    [[nodiscard]] size_type getNumberOfComponents()
        const noexcept override final;
    [[nodiscard]] bool isUniform() const noexcept override final;
    [[nodiscard]] std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx,
                    const real,
                    const real) const noexcept override final;
    [[nodiscard]] std::optional<PartialQuadratureFunctionEvaluatorResult>
    evaluate(Context&, const real, const real) const noexcept override final;
    //! \brief destructor
    ~UniformScalarPartialQuadratureFunctionEvaluatorBase() noexcept override;

   protected:
    /*!
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<real> getValue(
        Context&, const real, const real) const noexcept = 0;
  };  // end of UniformScalarPartialQuadratureFunctionEvaluatorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTIONEVALUATORBASE_HXX */
