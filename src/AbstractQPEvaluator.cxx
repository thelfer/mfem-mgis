/*!
 * \file   src/AbstractQPEvaluator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   07/03/2026
 */

#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/AbstractQPEvaluator.hxx"

namespace mfem_mgis {

  QPEvaluatorResult::QPEvaluatorResult(
      const ImmutablePartialQuadratureFunctionView& v) noexcept
      : ImmutablePartialQuadratureFunctionView(v) {}

  QPEvaluatorResult::QPEvaluatorResult(PartialQuadratureFunction&& v) noexcept
      : f(std::move(v)) {
    this->qspace = this->f->getPartialQuadratureSpacePointer();
    this->data_begin = this->f->getDataOffset();
    this->data_size = this->f->getNumberOfComponents();
    this->data_stride = this->f->getDataStride();
    this->immutable_values = this->f->getValues();
  }  // end of QPEvaluatorResult

  QPEvaluatorResult::QPEvaluatorResult(
      const PartialQuadratureFunction& v) noexcept
      : f(v) {
    this->qspace = this->f->getPartialQuadratureSpacePointer();
    this->data_begin = this->f->getDataOffset();
    this->data_size = this->f->getNumberOfComponents();
    this->data_stride = this->f->getDataStride();
    this->immutable_values = this->f->getValues();
  }  // end of QPEvaluatorResult

  QPEvaluatorResult::QPEvaluatorResult(QPEvaluatorResult&&) noexcept = default;

  QPEvaluatorResult::QPEvaluatorResult(const QPEvaluatorResult&) noexcept =
      default;

  QPEvaluatorResult& QPEvaluatorResult::operator=(
      QPEvaluatorResult&&) noexcept = default;

  QPEvaluatorResult::~QPEvaluatorResult() noexcept = default;

  AbstractQPEvaluator::~AbstractQPEvaluator() noexcept = default;

  std::optional<QPEvaluatorResult> evaluate(
      Context& ctx,
      const AbstractQPEvaluator& e,
      const real t,
      const real dt,
      const QPEvaluationOptions& opts) noexcept {
    auto oresult = e.evaluate(ctx, t, dt);
    if (isInvalid(oresult)) {
      return {};
    }
    if (!opts.ensure_contiguous_storage) {
      return oresult;
    }
    if (oresult->getNumberOfComponents() == oresult->getDataStride()) {
      // data are stored contigously in memory, nothing to be done
      if (oresult->getDataOffset() != 0) {
        return ctx.registerErrorMessage("invalid result");
      }
      return oresult;
    }
    // copy the result of the evaluation to a new partial quadrature function
    auto r =
        PartialQuadratureFunction(oresult->getPartialQuadratureSpacePointer(),
                                  oresult->getNumberOfComponents());
    r = *oresult;
    return {std::move(r)};
  }  // end of evaluate

}  // end of namespace mfem_mgis
