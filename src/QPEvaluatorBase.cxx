
/*!
 * \file   QPEvaluatorBase.cxx
 * \brief  This file implements the `QPEvaluatorBase`
 * class \author Thomas Helfer \date   12/03/2026
 */

#include "MFEMMGIS/QPEvaluatorBase.hxx"

namespace mfem_mgis {

  QPEvaluatorBase::QPEvaluatorBase(
      std::shared_ptr<const PartialQuadratureSpace> s)
      : qspace(std::move(s)) {
    if (this->qspace.get() == nullptr) {
      raise("invalid partial quadrature space");
    }
  }  // end of QPEvaluatorBase

  const PartialQuadratureSpace& QPEvaluatorBase::getQuadratureSpace()
      const noexcept {
    return *(this->qspace);
  }  // end of getQuadratureSpace

  std::shared_ptr<const PartialQuadratureSpace>
  QPEvaluatorBase::getPartialQuadratureSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getPartialQuadratureSpacePointer

  bool QPEvaluatorBase::isUniform() const noexcept {
    return false;
  }  // end of isUniform

  std::optional<std::variant<real, std::vector<real>>>
  QPEvaluatorBase::getUniformValue(Context& ctx,
                                   const real,
                                   const real) const noexcept {
    return ctx.registerErrorMessage("invalid call: evaluator is not valid");
  }  // end of getUniformValue

  QPEvaluatorBase::~QPEvaluatorBase() noexcept = default;

  size_type UniformScalarQPEvaluatorBase::getNumberOfComponents()
      const noexcept {
    return 1;
  }  // end of getNumberOfComponents

  bool UniformScalarQPEvaluatorBase::isUniform() const noexcept {
    return true;
  }  // end of isUniform()

  std::optional<std::variant<real, std::vector<real>>>
  UniformScalarQPEvaluatorBase::getUniformValue(Context& ctx,
                                                const real t,
                                                const real dt) const noexcept {
    return this->getValue(ctx, t, dt);
  }  // end of isUniform()

  std::optional<QPEvaluatorResult> UniformScalarQPEvaluatorBase::evaluate(
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

  UniformScalarQPEvaluatorBase::~UniformScalarQPEvaluatorBase() noexcept =
      default;

}  // end of namespace mfem_mgis
