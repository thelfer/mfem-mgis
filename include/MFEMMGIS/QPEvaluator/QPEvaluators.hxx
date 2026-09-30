/*!
 * \file   MFEMMGIS/QPEvaluator/QPEvaluators.hxx
 * \brief  This file declares a list of standard evaluators
 * \author Thomas HElfer
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

  //! \brief evaluator based on an uniform constant value
  struct MFEM_MGIS_EXPORT UniformConstantScalarQPEvaluator final
      : UniformScalarQPEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] r: value
     */
    UniformConstantScalarQPEvaluator(
        std::shared_ptr<const PartialQuadratureSpace> s, const real v);
    //! \brief destructor
    ~UniformConstantScalarQPEvaluator() noexcept override;

   protected:
    [[nodiscard]] std::optional<real> getValue(
        Context& ctx, const real t, const real dt) const noexcept override;

   private:
    const real value;
  };

  //! \brief evaluator based on an uniform value
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
    std::function<std::optional<real>(Context&, const real, const real)> fct;
  };

  //! \brief evaluator based on an uniform value
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
   * \return an evaluator of the gradient of the given name
   *
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the gradient variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator> makeGradientEvaluator(
      Context& ctx,
      const Material& m,
      std::string_view n,
      const TimeStepStage ts) noexcept;
  /*!
   * \return an evaluator of the thermodynamic force of the given name
   *
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator> makeThermodynamicForceEvaluator(
      Context& ctx,
      const Material& m,
      std::string_view n,
      const TimeStepStage ts) noexcept;
  /*!
   * \return an evaluator of the internal state variable of the given name
   *
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractQPEvaluator> makeInternalStateVariableEvaluator(
      Context& ctx,
      const Material& m,
      std::string_view n,
      const TimeStepStage ts) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_QPEVALUATOR_QPEVALUATORS_HXX */
