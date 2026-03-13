/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionEvaluators.hxx
 * \brief  This file declares a list of standard evaluators
 * \author Thomas HElfer
 * \date   12/03/2026
 */

#ifndef LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONEVALUATORS_HXX
#define LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONEVALUATORS_HXX

#include <memory>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionEvaluatorBase.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Material;

  //! \brief evaluator based on an uniform constant value
  struct MFEM_MGIS_EXPORT
      UniformConstantScalarPartialQuadratureFunctionEvaluator final
      : UniformScalarPartialQuadratureFunctionEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] r: value
     */
    UniformConstantScalarPartialQuadratureFunctionEvaluator(
        std::shared_ptr<const PartialQuadratureSpace>, const real);
    //! \brief destructor
    ~UniformConstantScalarPartialQuadratureFunctionEvaluator() noexcept
        override;

   protected:
    [[nodiscard]] std::optional<real> getValue(
        Context&, const real, const real) const noexcept override;

   private:
    const real value;
  };

  //! \brief evaluator based on an uniform value
  struct MFEM_MGIS_EXPORT UniformScalarPartialQuadratureFunctionEvaluator final
      : UniformScalarPartialQuadratureFunctionEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] f: function of time
     * \param[in] ts: time step stage
     */
    UniformScalarPartialQuadratureFunctionEvaluator(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::function<real(const real)>,
        const TimeStepStage);
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] f: function of time
     * \param[in] ts: time step stage
     */
    UniformScalarPartialQuadratureFunctionEvaluator(
        std::shared_ptr<const PartialQuadratureSpace>,
        std::function<std::optional<real>(Context&, const real)>,
        const TimeStepStage);
    //! \brief destructor
    ~UniformScalarPartialQuadratureFunctionEvaluator() noexcept override;

   protected:
    [[nodiscard]] std::optional<real> getValue(
        Context&, const real, const real) const noexcept override;

   private:
    std::function<std::optional<real>(Context&, const real, const real)> fct;
  };

  //! \brief evaluator based on an uniform value
  struct MFEM_MGIS_EXPORT StandardPartialQuadratureFunctionEvaluator final
      : PartialQuadratureFunctionEvaluatorBase {
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
    StandardPartialQuadratureFunctionEvaluator(
        std::shared_ptr<const PartialQuadratureSpace>,
        size_type,
        FirstFunctionType);
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] nc: number of components
     * \param[in] f: function
     */
    StandardPartialQuadratureFunctionEvaluator(
        std::shared_ptr<const PartialQuadratureSpace>,
        size_type,
        SecondFunctionType);
    //
    [[nodiscard]] size_type getNumberOfComponents()
        const noexcept override final;
    [[nodiscard]] std::optional<PartialQuadratureFunctionEvaluatorResult>
    evaluate(Context&, const real, const real) const noexcept override final;
    //! \brief destructor
    ~StandardPartialQuadratureFunctionEvaluator() noexcept override;

   private:
    //! \brief number of components
    const size_type n;
    //! \brief function
    std::variant<FirstFunctionType, SecondFunctionType> fct;
  };

  /*!
   * \return an evaluator of the internal state variable of the given name
   *
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  MFEM_MGIS_EXPORT [[nodiscard]]  //
  std::shared_ptr<AbstractPartialQuadratureFunctionEvaluator>
  makeInternalStateVariableEvaluator(Context&,
                                     const Material&,
                                     std::string_view,
                                     const TimeStepStage) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONEVALUATORS_HXX */
