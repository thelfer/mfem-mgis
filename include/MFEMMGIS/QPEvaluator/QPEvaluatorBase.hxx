/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluatorBase.hxx
 * \brief  This file declares the `QPEvaluatorBase` class
 * \author Thomas Helfer
 * \date   12/03/2026
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORBASE_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORBASE_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for most partial quadrature evaluators
   *
   * \note by default, evaluators are assumed to be spatially variable, thus the
   * default implementation of `isUniform` returns false and `getUniformValue`
   * returns an error.
   */
  struct MFEM_MGIS_EXPORT QPEvaluatorBase : AbstractQPEvaluator {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     */
    QPEvaluatorBase(std::shared_ptr<const PartialQuadratureSpace> s);
    //
    [[nodiscard]] const PartialQuadratureSpace& getQuadratureSpace()
        const noexcept override;
    [[nodiscard]] std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept override;
    //! \return if the evaluator is spatially uniform, always false
    [[nodiscard]] bool isUniform() const noexcept override;
    /*!
     * \brief report an error, the evaluator is not uniform
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step, unused
     * \param[in] dt: time increment, unused
     * \return an empty value
     */
    [[nodiscard]] std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx,
                    const real t,
                    const real dt) const noexcept override;
    //! \brief destructor
    ~QPEvaluatorBase() noexcept override;

   protected:
    //! \brief partial quadrature space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
  };  // end of struct QPEvaluatorBase

  //! \brief base class for evaluators returning uniform scalar values
  struct MFEM_MGIS_EXPORT UniformScalarQPEvaluatorBase : QPEvaluatorBase {
    using QPEvaluatorBase::QPEvaluatorBase;
    //
    //! \return the number of components, always 1
    [[nodiscard]] size_type getNumberOfComponents()
        const noexcept override final;
    //! \return if the evaluator is spatially uniform, always true
    [[nodiscard]] bool isUniform() const noexcept override final;
    [[nodiscard]] std::optional<std::variant<real, std::vector<real>>>
    getUniformValue(Context& ctx,
                    const real t,
                    const real dt) const noexcept override final;
    /*!
     * \brief evaluate the values at the integration points
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return a function equal to the uniform value at each integration point
     */
    [[nodiscard]] std::optional<QPEvaluatorResult> evaluate(
        Context& ctx,
        const real t,
        const real dt) const noexcept override final;
    //! \brief destructor
    ~UniformScalarQPEvaluatorBase() noexcept override;

   protected:
    /*!
     * \brief return the value of the evaluator
     * \return the value of the evaluator
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] virtual std::optional<real> getValue(
        Context& ctx, const real t, const real dt) const noexcept = 0;
  };  // end of UniformScalarQPEvaluatorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORBASE_HXX */
