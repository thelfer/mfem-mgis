/*!
 * \file   src/AbstractQPEvaluator.cxx
 * \brief  This file implements the classes and functions declared in
 * `MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx`
 * \author Thomas Helfer
 * \date   07/03/2026
 */

#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"

namespace mfem_mgis {

  QPEvaluatorResultStorage::QPEvaluatorResultStorage() noexcept = default;

  QPEvaluatorResultStorage::QPEvaluatorResultStorage(
      QPEvaluatorResultStorage&&) noexcept = default;

  QPEvaluatorResultStorage::QPEvaluatorResultStorage(
      PartialQuadratureFunction&& v) noexcept
      : f(std::move(v)) {}  // end of QPEvaluatorResultStorage

  QPEvaluatorResult QPEvaluatorResult::fromView(
      const ImmutablePartialQuadratureFunctionView& v) noexcept {
    return QPEvaluatorResult{v};
  }  // end of QPEvaluatorResult

  QPEvaluatorResult QPEvaluatorResult::fromFunction(
      PartialQuadratureFunction&& v) noexcept {
    return QPEvaluatorResult{std::move(v)};
  }

  QPEvaluatorResult::QPEvaluatorResult(
      const ImmutablePartialQuadratureFunctionView& v) noexcept
      : ImmutablePartialQuadratureFunctionView(v) {}

  QPEvaluatorResult::QPEvaluatorResult(PartialQuadratureFunction&& v) noexcept
      : QPEvaluatorResultStorage(std::move(v)),
        ImmutablePartialQuadratureFunctionView(*(this->f)) {
  }  // end of QPEvaluatorResult

  QPEvaluatorResult::QPEvaluatorResult(QPEvaluatorResult&&) noexcept = default;

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
    auto of = PartialQuadratureFunction::copy(ctx, *oresult);
    if (isInvalid(of)) {
      return {};
    }
    return QPEvaluatorResult::fromFunction(std::move(*of));
  }  // end of evaluate

}  // end of namespace mfem_mgis
