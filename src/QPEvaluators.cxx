/*!
 * \file   QPEvaluators.cxx
 * \brief  This file implements some of the most standard partial quadrature
 * function evaluators
 * \author Thomas Helfer
 * \date   12/03/2026
 */

#include <algorithm>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/QPEvaluators.hxx"

namespace mfem_mgis {

  UniformConstantScalarQPEvaluator::UniformConstantScalarQPEvaluator(
      std::shared_ptr<const PartialQuadratureSpace> s, const real v)
      : UniformScalarQPEvaluatorBase(s),
        value(v) {}  // end of UniformConstantScalarQPEvaluator

  std::optional<real> UniformConstantScalarQPEvaluator::getValue(
      Context&, const real, const real) const noexcept {
    return this->value;
  }  // end of UniformConstantScalarQPEvaluator

  UniformConstantScalarQPEvaluator::
      ~UniformConstantScalarQPEvaluator() noexcept = default;

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

  UniformScalarQPEvaluator::UniformScalarQPEvaluator(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<real(const real)> f,
      const TimeStepStage ts)
      : UniformScalarQPEvaluatorBase(s),
        fct(buildFunction(throwing, f, ts)) {
  }  // end of UniformScalarQPEvaluator

  UniformScalarQPEvaluator::UniformScalarQPEvaluator(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::function<std::optional<real>(Context&, const real)> f,
      const TimeStepStage ts)
      : UniformScalarQPEvaluatorBase(s),
        fct(buildFunction(throwing, f, ts)) {
  }  // end of UniformScalarQPEvaluator

  std::optional<real> UniformScalarQPEvaluator::getValue(
      Context& ctx, const real t, const real dt) const noexcept {
    return this->fct(ctx, t, dt);
  }  // end of UniformScalarQPEvaluator

  UniformScalarQPEvaluator::~UniformScalarQPEvaluator() noexcept = default;

  StandardQPEvaluator::StandardQPEvaluator(
      std::shared_ptr<const PartialQuadratureSpace> s,
      size_type nc,
      FirstFunctionType f)
      : QPEvaluatorBase(s), n(nc), fct(f) {
    if (!f) {
      raise("invalid function");
    }
  }  // end of StandardQPEvaluator

  StandardQPEvaluator::StandardQPEvaluator(
      std::shared_ptr<const PartialQuadratureSpace> s,
      size_type nc,
      SecondFunctionType f)
      : QPEvaluatorBase(s), n(nc), fct(f) {
    if (!f) {
      raise("invalid function");
    }
  }  // end of StandardQPEvaluator

  size_type StandardQPEvaluator::getNumberOfComponents() const noexcept {
    return this->n;
  }  // end of StandardQPEvaluator

  std::optional<QPEvaluatorResult> StandardQPEvaluator::evaluate(
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

  StandardQPEvaluator::~StandardQPEvaluator() noexcept = default;

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

  std::shared_ptr<AbstractQPEvaluator> makeGradientEvaluator(
      Context& ctx,
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
    return make_shared<StandardQPEvaluator>(ctx, qspace, *os, extract);
  }  // end of makeGradientEvaluator

  std::shared_ptr<AbstractQPEvaluator> makeThermodynamicForceEvaluator(
      Context& ctx,
      const Material& m,
      std::string_view n,
      const TimeStepStage ts) noexcept {
    using namespace mgis::behaviour;
    const auto qspace = m.getPartialQuadratureSpacePointer();
    const auto os =
        getVariableSize(ctx, m.b.thermodynamic_forces, n, m.b.hypothesis);
    if (isInvalid(os)) {
      return {};
    }
    // make a copy of the name, so that the lambda also take a copy
    auto name = std::string{n};
    auto extract = [&m, ts, name](Context& ectx, const real,
                                  const real) noexcept {
      return getThermodynamicForce(ectx, m, name, ts);
    };
    return make_shared<StandardQPEvaluator>(ctx, qspace, *os, extract);
  }  // end of makeThermodynamicForceEvaluator

  std::shared_ptr<AbstractQPEvaluator> makeInternalStateVariableEvaluator(
      Context& ctx,
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
    return make_shared<StandardQPEvaluator>(ctx, qspace, *os, extract);
  }  // end of makeInternalStateVariableEvaluator

}  // end of namespace mfem_mgis
