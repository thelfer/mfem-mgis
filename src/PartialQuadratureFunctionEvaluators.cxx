/*!
 * \file   PartialQuadratureFunctionEvaluators.cxx
 * \brief  This file implements some of the most standard partial quadrature
 * function evaluators
 * \author Thomas Helfer
 * \date   12/03/2026
 */

#include <algorithm>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionEvaluators.hxx"

namespace mfem_mgis {

  UniformConstantScalarPartialQuadratureFunctionEvaluator::
      UniformConstantScalarPartialQuadratureFunctionEvaluator(
          std::shared_ptr<const PartialQuadratureSpace> s, const real v)
      : UniformScalarPartialQuadratureFunctionEvaluatorBase(s),
        value(v) {
  }  // end of UniformConstantScalarPartialQuadratureFunctionEvaluator

  std::optional<real>
  UniformConstantScalarPartialQuadratureFunctionEvaluator::getValue(
      Context&, const real, const real) const noexcept {
    return this->value;
  }  // end of UniformConstantScalarPartialQuadratureFunctionEvaluator

  UniformConstantScalarPartialQuadratureFunctionEvaluator::
      ~UniformConstantScalarPartialQuadratureFunctionEvaluator() noexcept =
          default;

  static std::function<std::optional<real>(Context&, const real, const real)>
  buildFunction(attributes::Throwing,
                std::function<real(const real)> f,
                const TimeStepStage ts) {
    if (!f) {
      raise("invalid function");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      return
          [f](Context& ctx, const real t, const real) -> std::optional<real> {
            try {
              return f(t);
            } catch (...) {
              std::ignore = registerExceptionInErrorBacktrace(ctx);
            }
            return {};
          };
    }
    return
        [f](Context& ctx, const real t, const real dt) -> std::optional<real> {
          try {
            return f(t + dt);
          } catch (...) {
            std::ignore = registerExceptionInErrorBacktrace(ctx);
          }
          return {};
        };
  }  // end of buildFunction

  static std::function<std::optional<real>(Context&, const real, const real)>
  buildFunction(attributes::Throwing,
                std::function<std::optional<real>(Context& ctx, const real)> f,
                const TimeStepStage ts) {
    if (!f) {
      raise("invalid function");
    }
    if (ts == TimeStepStage::BEGINNING_OF_TIME_STEP) {
      return
          [f](Context& ctx, const real t, const real) -> std::optional<real> {
            try {
              return f(ctx, t);
            } catch (...) {
              std::ignore = registerExceptionInErrorBacktrace(ctx);
            }
            return {};
          };
    }
    return
        [f](Context& ctx, const real t, const real dt) -> std::optional<real> {
          try {
            return f(ctx, t + dt);
          } catch (...) {
            std::ignore = registerExceptionInErrorBacktrace(ctx);
          }
          return {};
        };
  }  // end of buildFunction

  UniformScalarPartialQuadratureFunctionEvaluator::
      UniformScalarPartialQuadratureFunctionEvaluator(
          std::shared_ptr<const PartialQuadratureSpace> s,
          std::function<real(const real)> f,
          const TimeStepStage ts)
      : UniformScalarPartialQuadratureFunctionEvaluatorBase(s),
        fct(buildFunction(throwing, f, ts)) {
  }  // end of UniformScalarPartialQuadratureFunctionEvaluator

  UniformScalarPartialQuadratureFunctionEvaluator::
      UniformScalarPartialQuadratureFunctionEvaluator(
          std::shared_ptr<const PartialQuadratureSpace> s,
          std::function<std::optional<real>(Context&, const real)> f,
          const TimeStepStage ts)
      : UniformScalarPartialQuadratureFunctionEvaluatorBase(s),
        fct(buildFunction(throwing, f, ts)) {
  }  // end of UniformScalarPartialQuadratureFunctionEvaluator

  std::optional<real> UniformScalarPartialQuadratureFunctionEvaluator::getValue(
      Context& ctx, const real t, const real dt) const noexcept {
    return this->fct(ctx, t, dt);
  }  // end of UniformScalarPartialQuadratureFunctionEvaluator

  UniformScalarPartialQuadratureFunctionEvaluator::
      ~UniformScalarPartialQuadratureFunctionEvaluator() noexcept = default;

  StandardPartialQuadratureFunctionEvaluator::
      StandardPartialQuadratureFunctionEvaluator(
          std::shared_ptr<const PartialQuadratureSpace> s,
          size_type nc,
          FirstFunctionType f)
      : PartialQuadratureFunctionEvaluatorBase(s), n(nc), fct(f) {
    if (!f) {
      raise("invalid function");
    }
  }  // end of StandardPartialQuadratureFunctionEvaluator

  StandardPartialQuadratureFunctionEvaluator::
      StandardPartialQuadratureFunctionEvaluator(
          std::shared_ptr<const PartialQuadratureSpace> s,
          size_type nc,
          SecondFunctionType f)
      : PartialQuadratureFunctionEvaluatorBase(s), n(nc), fct(f) {
    if (!f) {
      raise("invalid function");
    }
  }  // end of StandardPartialQuadratureFunctionEvaluator

  size_type StandardPartialQuadratureFunctionEvaluator::getNumberOfComponents()
      const noexcept {
    return this->n;
  }  // end of StandardPartialQuadratureFunctionEvaluator

  std::optional<PartialQuadratureFunctionEvaluatorResult>
  StandardPartialQuadratureFunctionEvaluator::evaluate(
      Context& ctx, const real t, const real dt) const noexcept {
    try {
      if (std::holds_alternative<FirstFunctionType>(this->fct)) {
        const auto f = std::get<FirstFunctionType>(this->fct);
        const auto ovalue = f(ctx, t, dt);
        if (isInvalid(ovalue)) {
          return {};
        }
        if (ovalue->getNumberOfComponents() != this->n) {
          return ctx.registerErrorMessage(
              "invalid number of components (" +
              std::to_string(ovalue->getNumberOfComponents()) + ", expected " +
              std::to_string(this->n) + ")");
        }
        return *ovalue;
      }
      const auto f = std::get<SecondFunctionType>(this->fct);
      const auto ovalue = f(ctx, t, dt);
      if (isInvalid(ovalue)) {
        return {};
      }
      if (ovalue->getNumberOfComponents() != this->n) {
        return ctx.registerErrorMessage(
            "invalid number of components (" +
            std::to_string(ovalue->getNumberOfComponents()) + ", expected " +
            std::to_string(this->n) + ")");
      }
      return *ovalue;
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of evaluate

  StandardPartialQuadratureFunctionEvaluator::
      ~StandardPartialQuadratureFunctionEvaluator() noexcept = default;

  [[nodiscard]] static std::optional<size_type> getVariableSize(
      Context& ctx,
      const std::vector<mgis::behaviour::Variable>& variables,
      std::string_view n,
      const mgis::behaviour::Hypothesis h) {
    if (!contains(variables, n)) {
      return ctx.registerErrorMessage("no variable named '" + std::string{n} +
                                      "'");
    }
    const auto ov = getVariable(ctx, variables, n);
    if (isInvalid(ov)) {
      return {};
    }
    return getVariableSize(ctx, *ov, h);
  }  // end of getVariableSize

  std::shared_ptr<AbstractPartialQuadratureFunctionEvaluator>
  makeGradientEvaluator(Context& ctx,
                        const Material& m,
                        std::string_view n,
                        const TimeStepStage ts) noexcept {
    using namespace mgis::behaviour;
    const auto qspace = m.getPartialQuadratureSpacePointer();
    const auto os = getVariableSize(ctx, m.b.gradients, n, m.b.hypothesis);
    if (isInvalid(os)) {
      return {};
    }
    // make a copy of the name, so that the lambda also take a copy
    auto name = std::string{n};
    auto extract = [&m, ts, name](Context& ectx, const real,
                                  const real) noexcept {
      return getGradient(ectx, m, name, ts);
    };
    return make_shared<StandardPartialQuadratureFunctionEvaluator>(
        ctx, qspace, *os, extract);
  }  // end of makeGradientEvaluator

  std::shared_ptr<AbstractPartialQuadratureFunctionEvaluator>
  makeInternalStateVariableEvaluator(Context& ctx,
                                     const Material& m,
                                     std::string_view n,
                                     const TimeStepStage ts) noexcept {
    using namespace mgis::behaviour;
    const auto qspace = m.getPartialQuadratureSpacePointer();
    const auto os = getVariableSize(ctx, m.b.isvs, n, m.b.hypothesis);
    if (isInvalid(os)) {
      return {};
    }
    // make a copy of the name, so that the lambda also take a copy
    auto name = std::string{n};
    auto extract = [&m, ts, name](Context& ectx, const real,
                                  const real) noexcept {
      return getInternalStateVariable(ectx, m, name, ts);
    };
    return make_shared<StandardPartialQuadratureFunctionEvaluator>(
        ctx, qspace, *os, extract);
  }  // end of makeInternalStateVariableEvaluator

}  // end of namespace mfem_mgis