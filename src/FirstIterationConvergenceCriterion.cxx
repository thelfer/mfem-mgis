/*!
 * \file   src/FirstIterationConvergenceCriterion.cxx
 * \brief  This file implements the FirstIterationConvergenceCriterion class
 * \date   02/09/2024
 */

#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FirstIterationConvergenceCriterion.hxx"

namespace mfem_mgis {

  std::string FirstIterationConvergenceCriterion::getDescription() noexcept {
    return "a convergence criterion which checks if every coupling item "
           "converges at the first iteration";
  }  // end of getDescription

  std::map<std::string, std::string>
  FirstIterationConvergenceCriterion::getParametersDescription() noexcept {
    return {};
  }  // end of getParametersDescription

  FirstIterationConvergenceCriterion::FirstIterationConvergenceCriterion() =
      default;

  bool FirstIterationConvergenceCriterion::
      performInitializationTaksAtTheBeginningOfTheTimeStep(
          Context &, const TimeStep &) noexcept {
    return true;
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  [[nodiscard]] static std::optional<bool> treatStandardCouplingItem(
      Context &ctx, const Parameters &itemOutput) noexcept {
    const auto oname = get<std::string>(ctx, itemOutput, "Name");
    const auto oout = get<Parameters>(ctx, itemOutput, "Output");
    if (!(areValid(oname, oout))) {
      return {};
    }
    if (contains(*oout, "Solver")) {
      // models based on a solver
      const auto &osolver = get<Parameters>(ctx, *oout, "Solver");
      if (isInvalid(osolver)) {
        return {};
      }
      if (!contains(*osolver, "NumberOfIterations")) {
        return true;
      }
      const auto on = get<size_type>(ctx, *osolver, "NumberOfIterations");
      if (isInvalid(on)) {
        return {};
      }
      return *on == 0;
    }
    return true;
  }  // end of treatStandardCouplingItem

  static std::optional<bool> treatCouplingScheme(
      Context &ctx, const Parameters &out) noexcept {
    if (!contains(out, "ItemsOutputs")) {
      return true;
    }
    const auto ooutputs = get<std::vector<Parameter>>(ctx, out, "ItemsOutputs");
    if (isInvalid(ooutputs)) {
      return {};
    }
    auto success = true;
    for (const auto &o : *ooutputs) {
      const auto r = [&ctx, &o]() -> std::optional<bool> {
        const auto &oitem_output = get<Parameters>(ctx, o);
        if (!oitem_output) {
          return {};
        }
        if (contains(*oitem_output, "ItemsOutputs")) {
          // coupling scheme
          return treatCouplingScheme(ctx, *oitem_output);
        }
        return treatStandardCouplingItem(ctx, *oitem_output);
      }();
      if (isInvalid(r)) {
        return {};
      }
      success = success && *r;
    }
    return success;
  }

  std::optional<bool> FirstIterationConvergenceCriterion::check(
      Context &ctx, const ComputeNextStateOutput &output) const noexcept {
    return treatCouplingScheme(ctx, output);
  }  // end of check

  bool FirstIterationConvergenceCriterion::update(Context &) noexcept {
    return true;
  }

  bool FirstIterationConvergenceCriterion::revert(Context &) noexcept {
    return true;
  }

  FirstIterationConvergenceCriterion::
      ~FirstIterationConvergenceCriterion() noexcept = default;

}  // namespace mfem_mgis