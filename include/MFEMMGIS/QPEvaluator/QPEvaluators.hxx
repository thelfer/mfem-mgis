/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluators.hxx
 * \brief  This file declares a list of standard evaluators
 * \author Thomas Helfer
 * \date   12/03/2026
 */

#ifndef LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORS_HXX
#define LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORS_HXX

#include <memory>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/QPEvaluator/QPEvaluatorBase.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Material;

  //! \brief evaluator based on a uniform constant value
  struct MFEM_MGIS_EXPORT UniformConstantScalarQPEvaluator final
      : UniformScalarQPEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] v: value
     */
    UniformConstantScalarQPEvaluator(
        std::shared_ptr<const PartialQuadratureSpace> s, const real v);
    //! \brief destructor
    ~UniformConstantScalarQPEvaluator() noexcept override;

   protected:
    /*!
     * \brief return the constant value
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step, unused
     * \param[in] dt: time increment, unused
     * \return the value
     */
    [[nodiscard]] std::optional<real> getValue(
        Context& ctx, const real t, const real dt) const noexcept override;

   private:
    //! \brief uniform value
    const real value;
  };

  //! \brief evaluator based on a uniform value
  struct MFEM_MGIS_EXPORT UniformScalarQPEvaluator final
      : UniformScalarQPEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] f: function of time
     * \param[in] ts: time step stage
     */
    UniformScalarQPEvaluator(std::shared_ptr<const PartialQuadratureSpace> s,
                             std::function<real(const real)> f,
                             const TimeStepStage ts);
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] f: function of time
     * \param[in] ts: time step stage
     */
    UniformScalarQPEvaluator(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::function<std::optional<real>(Context&, const real)> f,
        const TimeStepStage ts);
    //! \brief destructor
    ~UniformScalarQPEvaluator() noexcept override;

   protected:
    [[nodiscard]] std::optional<real> getValue(
        Context& ctx, const real t, const real dt) const noexcept override;

   private:
    //! \brief function returning the uniform value
    std::function<std::optional<real>(Context&, const real, const real)> fct;
  };

  //! \brief evaluator based on a function returning the values
  struct MFEM_MGIS_EXPORT StandardQPEvaluator final : QPEvaluatorBase {
    //! \brief a simple alias
    using FirstFunctionType =
        std::function<std::optional<ImmutablePartialQuadratureFunctionView>(
            Context&, const real, const real)>;
    //! \brief a simple alias
    using SecondFunctionType =
        std::function<std::optional<PartialQuadratureFunction>(
            Context&, const real, const real)>;
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] nc: number of components
     * \param[in] f: function
     */
    StandardQPEvaluator(std::shared_ptr<const PartialQuadratureSpace> s,
                        size_type nc,
                        FirstFunctionType f);
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] nc: number of components
     * \param[in] f: function
     */
    StandardQPEvaluator(std::shared_ptr<const PartialQuadratureSpace> s,
                        size_type nc,
                        SecondFunctionType f);
    //
    [[nodiscard]] size_type getNumberOfComponents()
        const noexcept override final;
    /*!
     * \brief evaluate the values at the integration points
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return the values returned by the function
     */
    [[nodiscard]] std::optional<QPEvaluatorResult> evaluate(
        Context& ctx,
        const real t,
        const real dt) const noexcept override final;
    //! \brief destructor
    ~StandardQPEvaluator() noexcept override;

   private:
    //! \brief number of components
    const size_type n;
    //! \brief function
    std::variant<FirstFunctionType, SecondFunctionType> fct;
  };

  /*!
   * \brief create an evaluator of a gradient
   * \return an evaluator of the gradient of the given name
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the gradient variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator>
  makeGradientEvaluator(Context& ctx,
                        const Material& m,
                        std::string_view n,
                        const TimeStepStage ts) noexcept;
  /*!
   * \brief create an evaluator of a thermodynamic force
   * \return an evaluator of the thermodynamic force of the given name
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator>
  makeThermodynamicForceEvaluator(Context& ctx,
                                  const Material& m,
                                  std::string_view n,
                                  const TimeStepStage ts) noexcept;
  /*!
   * \brief create an evaluator of an internal state variable
   * \return an evaluator of the internal state variable of the given name
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator>
  makeInternalStateVariableEvaluator(Context& ctx,
                                     const Material& m,
                                     std::string_view n,
                                     const TimeStepStage ts) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORS_HXX */
