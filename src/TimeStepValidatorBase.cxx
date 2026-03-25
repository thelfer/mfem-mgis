
/*!
 * \file   src/TimeStepValidatorBase.cxx
 * \brief  This class implements the `TimeStepValidatorBase` class
 * \date   04/12/2023
 */

#include <limits>
#include <iostream>
#include <algorithm>
#include "MFEMMGIS/TimeStepValidatorBase.hxx"

namespace mfem_mgis {

  TimeStepValidatorBase::TimeStepValidatorBase() noexcept = default;

  void TimeStepValidatorBase::addValidator(
      const ExternalValidator &v) noexcept {
    this->addValidator(
        "validator " + std::to_string(this->externalValidators.size()), v);
  }  // end of addValidator

  void TimeStepValidatorBase::addValidator(
      std::string_view n, const ExternalValidator &v) noexcept {
    this->externalValidators.push_back({std::string{n}, v});
  }  // end of addValidator

  std::optional<AbstractTimeStepValidator::Result>
  TimeStepValidatorBase::callExternalValidators(Context &ctx) const noexcept {
    auto r = Result{};
    for (const auto &[n, v] : this->externalValidators) {
      const auto r2 = [&ctx, &v]() -> std::optional<std::pair<bool, double>> {
        try {
          return std::invoke(v);
        } catch (...) {
          std::ignore = mgis::registerExceptionInErrorBacktrace(ctx);
        }
        return {};
      }();
      if (!r2.has_value()) {
        return ctx.registerErrorMessage(
            "calling external time step validator '" + n + "' failed");
      }
      r.isValid = r.isValid && r2->first;
      r.recommendedTimeIncrement =
          std::min(r.recommendedTimeIncrement, r2->second);
      if (!r2->first) {
        r.reasons.push_back(
            "time step rejected by external time step validator '" + n + "'");
      }
    }
    return r;
  }  // end of callExternalValidators_

  TimeStepValidatorBase::~TimeStepValidatorBase() = default;

}  // namespace mfem_mgis
