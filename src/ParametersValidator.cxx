/*!
 * \file   src/ParametersValidator.cxx
 * \brief  This files implements the `ParametersValidator` class
 * \author Thomas Helfer
 * \date   19/09/2026
 */

#include <algorithm>
#include "MFEMMGIS/Utilities/ParametersValidator.hxx"

namespace mfem_mgis {

  InvalidResult ParametersValidator::reportUnmatchedTypeError(
      Context& ctx, const std::string& k) noexcept {
    return ctx.registerErrorMessage("parameter '" + k +
                                    "' does not hold the expected type");
  }  // end of reportUnmatchedTypeError

  InvalidResult ParametersValidator::reportUnmatchedTypesError(
      Context& ctx, const std::string& k) noexcept {
    return ctx.registerErrorMessage(
        "parameter '" + k + "' does not hold any of the expected types");
  }  // end of reportUnmatchedTypesError

  ParametersValidator::ParametersValidator() noexcept = default;
  ParametersValidator::ParametersValidator(ParametersValidator&&) noexcept =
      default;
  ParametersValidator::ParametersValidator(
      const ParametersValidator&) noexcept = default;
  ParametersValidator& ParametersValidator::operator=(
      ParametersValidator&&) noexcept = default;
  ParametersValidator& ParametersValidator::operator=(
      const ParametersValidator&) noexcept = default;

  ParametersValidator& ParametersValidator::addStrictlyPositiveIntegerCheck(
      const std::string& k, const AddArguments& opts) noexcept {
    const auto v = [k](Context& ctx, const Parameter& d) noexcept -> bool {
      if (!is<int>(d)) {
        return ctx.registerErrorMessage("parameter '" + k +
                                        "' is not an integer");
      }
      if (get<int>(throwing, d) <= 0) {  // this can't throw
        return ctx.registerErrorMessage("parameter '" + k +
                                        "' is not strictly postive");
      }
      return true;
    };
    return this->add(k, v, opts);
  }  // end of addStrictlyPositiveIntegerCheck

  ParametersValidator& ParametersValidator::addStrictlyPositiveIntegerCheck(
      const std::string& k,
      const std::string& d,
      const AddArguments& opts) noexcept {
    this->add(k, d, opts);  // pre-declare the key with the documentation
    return this->addStrictlyPositiveIntegerCheck(k, opts);
  }  // end of addStrictlyPositiveIntegerCheck

  ParametersValidator& ParametersValidator::add(
      const std::vector<std::string>& keys, const AddArguments& opts) noexcept {
    for (const auto& k : keys) {
      this->add(k, opts);
    }
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::add(
      const std::string& k, const AddArguments& opts) noexcept {
    this->addKey(k, opts);
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::add(
      const std::string& k,
      const std::string& d,
      const AddArguments& opts) noexcept {
    this->addKey(k, opts);
    const auto p = this->allowed_keys.find(k);
    if (p != this->allowed_keys.end()) {
      // this branch is always executed
      if (p->second.empty()) {
        this->allowed_keys[k] = d;
      }
    }
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::add(
      const std::map<std::string, std::string>& m,
      const AddArguments& opts) noexcept {
    for (const auto& [k, d] : m) {
      this->add(k, d, opts);
    }
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::add(
      const std::string& k,
      const ParameterValidator& f,
      const AddArguments& opts) noexcept {
    this->addKey(k, opts);
    if (f) {  // ignoring invalid validator
      this->validators[k].push_back(f);
    }
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::add(
      const std::string& k,
      const std::string& d,
      const ParameterValidator& f,
      const AddArguments& opts) noexcept {
    this->add(k, d, opts);
    if (f) {  // ignoring invalid validator
      this->validators[k].push_back(f);
    }
    return *this;
  }  // end of add

  ParametersValidator& ParametersValidator::addIncompatibleParametersList(
      const std::vector<std::string>& keys, const AddArguments& opts) noexcept {
    for (const auto& k : keys) {
      this->addKey(k, {});
    }
    if (keys.size() > 1) {
      this->incompatibilities.push_back(keys);
    }
    if (opts.required) {
      if (keys.size() > 1) {
        this->required_keys_in_set.push_back(keys);
      } else if (keys.size() == 1) {
        this->required_keys.insert(keys.at(0));
      }
    }
    return *this;
  }  // end of addIncompatibleParametersList

  ParametersValidator& ParametersValidator::addIncompatibleParametersList(
      const std::map<std::string, std::string>& m,
      const AddArguments& opts) noexcept {
    auto keys = std::vector<std::string>{};
    for (const auto& [k, d] : m) {
      keys.push_back(k);
      this->add(k, d, {.required = false});
    }
    return this->addIncompatibleParametersList(keys, opts);
  }  // end of addIncompatibleParametersList

  void ParametersValidator::addKey(const std::string& k,
                                   const AddArguments& opts) noexcept {
    const auto p = this->allowed_keys.find(k);
    if (p == this->allowed_keys.end()) {
      this->allowed_keys.insert({k, ""});
    }
    if (opts.required) {
      this->required_keys.insert(k);
    }
  }  // end of addKey

  void ParametersValidator::validate(attributes::Throwing,
                                     const Parameters& m) const {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    this->validate(ctx, m) | or_raise;
  }  // end of validate

  bool ParametersValidator::validate(Context& ctx,
                                     const Parameters& m) const noexcept {
    for (const auto& keys : this->incompatibilities) {
      for (const auto& k : keys) {
        ctx.assertOrTerminate(
            !this->required_keys.contains(k),
            "parameter '" + k +
                "' is declared as being incompatible "
                "with some other parameters and is also declared "
                "required. This does not make sense");
      }
    }
    for (const auto& k : this->required_keys) {
      if (!m.contains(k)) {
        return ctx.registerErrorMessage("required parameter '" + k +
                                        "' is missing");
      }
    }
    for (const auto& keys : this->required_keys_in_set) {
      const auto found = [&keys, &m] {
        for (const auto& k : keys) {
          if (m.contains(k)) {
            return true;
          }
        }
        return false;
      }();
      if (!found) {
        auto msg =
            std::string{"one of the following parameter must be defined:"};
        auto first = true;
        for (const auto& k : keys) {
          if (first) {
            msg += " '" + k + '\'';
          } else {
            msg += ", '" + k + '\'';
          }
          first = false;
        }
        return ctx.registerErrorMessage(msg);
      }
    }
    for (const auto& [k, v] : m) {
      static_cast<void>(v);
      const auto p = this->allowed_keys.find(k);
      if (p == this->allowed_keys.end()) {
        auto msg = std::string("invalid parameter '" + k + "'. ");
        if (!this->allowed_keys.empty()) {
          msg += "Valid parameters are:";
          for (const auto& [k2, d2] : this->allowed_keys) {
            msg += "\n- '" + k2 + "'";
            if (!d2.empty()) {
              msg += ": " + d2;
            } else {
              msg += " (undocumented)";
            }
          }
        }
        return ctx.registerErrorMessage(msg);
      }
      for (const auto& keys : this->incompatibilities) {
        if (std::find(keys.begin(), keys.end(), k) != keys.end()) {
          for (const auto& k2 : keys) {
            if (k2 == k) {
              continue;
            }
            if (m.contains(k2)) {
              return ctx.registerErrorMessage(
                  "parameters '" + k + "' and '" + k2 +
                  +"' are exclusive: only one shall be defined");
            }
          }
        }
      }
      const auto pvs = this->validators.find(k);
      if (pvs == this->validators.end()) {
        continue;
      }
      for (const auto& validator : pvs->second) {
        try {
          if (!validator(ctx, v)) {
            return false;
          }
        } catch (std::exception& e) {
          return ctx.registerErrorMessage(
              "ParametersValidator::validate: "
              "invalid value for key '" +
              k + "' (" + std::string(e.what()) + ")");
        } catch (...) {
          return ctx.registerErrorMessage(
              "ParametersValidator::validate: "
              "invalid value for key '" +
              k + "' (unhandled exception was thrown)");
        }
      }
    }
    return true;
  }

  const std::map<std::string, std::string, std::less<>>&
  ParametersValidator::getAllowedParameters() const noexcept {
    return this->allowed_keys;
  }  // end of getAllowedParameters

  std::optional<std::string> ParametersValidator::getDescription(
      Context& ctx, std::string_view k) const noexcept {
    const auto p = this->allowed_keys.find(k);
    if (p != this->allowed_keys.end()) {
      return p->second;
    }
    return ctx.registerErrorMessage("parameter '" + std::string{k} +
                                    "' is not declared");
  }  // end of getDescription

  ParametersValidator::~ParametersValidator() = default;

}  // end of namespace mfem_mgis