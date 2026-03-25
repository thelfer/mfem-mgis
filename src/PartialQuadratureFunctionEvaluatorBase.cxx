
/*!
 * \file   PartialQuadratureFunctionEvaluatorBase.cxx
 * \brief  This file implements the `PartialQuadratureFunctionEvaluatorBase`
 * class \author Thomas Helfer \date   12/03/2026
 */

#include "MFEMMGIS/PartialQuadratureFunctionEvaluatorBase.hxx"

namespace mfem_mgis {

  PartialQuadratureFunctionEvaluatorBase::
      PartialQuadratureFunctionEvaluatorBase(
          std::shared_ptr<const PartialQuadratureSpace> s)
      : qspace(std::move(s)) {
    if (this->qspace.get() == nullptr) {
      raise("invalid partial quadrature space");
    }
  }  // end of PartialQuadratureFunctionEvaluatorBase

  const PartialQuadratureSpace&
  PartialQuadratureFunctionEvaluatorBase::getQuadratureSpace() const noexcept {
    return *(this->qspace);
  }  // end of getQuadratureSpace

  std::shared_ptr<const PartialQuadratureSpace>
  PartialQuadratureFunctionEvaluatorBase::getPartialQuadratureSpacePointer()
      const noexcept {
    return this->qspace;
  }  // end of getPartialQuadratureSpacePointer

  bool PartialQuadratureFunctionEvaluatorBase::isUniform() const noexcept {
    return false;
  }  // end of isUniform

  std::optional<std::variant<real, std::vector<real>>>
  PartialQuadratureFunctionEvaluatorBase::getUniformValue(
      Context& ctx, const real, const real) const noexcept {
    return ctx.registerErrorMessage("invalid call: evaluator is not valid");
  }  // end of getUniformValue

  PartialQuadratureFunctionEvaluatorBase::
      ~PartialQuadratureFunctionEvaluatorBase() noexcept = default;

  size_type
  UniformScalarPartialQuadratureFunctionEvaluatorBase::getNumberOfComponents()
      const noexcept {
    return 1;
  }  // end of getNumberOfComponents

  bool UniformScalarPartialQuadratureFunctionEvaluatorBase::isUniform()
      const noexcept {
    return true;
  }  // end of isUniform()

  std::optional<std::variant<real, std::vector<real>>>
  UniformScalarPartialQuadratureFunctionEvaluatorBase::getUniformValue(
      Context& ctx, const real t, const real dt) const noexcept {
    return this->getValue(ctx, t, dt);
  }  // end of isUniform()

  std::optional<PartialQuadratureFunctionEvaluatorResult>
  UniformScalarPartialQuadratureFunctionEvaluatorBase::evaluate(
      Context& ctx, const real t, const real dt) const noexcept {
    const auto ov = this->getValue(ctx, t, dt);
    if (isInvalid(ov)) {
      return {};
    }
    auto r = PartialQuadratureFunction{this->qspace};
    auto values = r.getValues();

    std::fill(values.begin(), values.end(), *ov);
    return {std::move(r)};
  }

  UniformScalarPartialQuadratureFunctionEvaluatorBase::
      ~UniformScalarPartialQuadratureFunctionEvaluatorBase() noexcept = default;

}  // end of namespace mfem_mgis
