/*!
 * \file   src/Parameter.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/03/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameter.hxx"

namespace mfem_mgis::internals {

  template <typename ValueType>
  [[nodiscard]] static std::optional<std::vector<ValueType>>
  convertToVector_impl(Context& ctx, const Parameter& p) noexcept {
    if (is<ValueType>(p)) {
      const auto& value = get<ValueType>(throwing, p);
      return std::vector<ValueType>(1, value);
    }
    if (!is<std::vector<Parameter>>(p)) {
      return ctx.registerErrorMessage(
          "can't convert parameter to vector: invalid parameter type");
    }
    const auto& values = get<std::vector<Parameter>>(throwing, p);
    auto r = std::vector<ValueType>{};
    r.reserve(values.size());
    for (const auto& v : values) {
      if (!is<ValueType>(v)) {
        return ctx.registerErrorMessage(
            "can't convert parameter to vector: one of the element is not "
            "convertible to the expected type");
      }
      r.push_back(get<ValueType>(throwing, v));
    }
    return r;
  }  // end of convertToVector_impl

  template <>
  std::optional<std::vector<real>> convertToVector<real>(
      Context& ctx, const Parameter& p) noexcept {
    return convertToVector_impl<real>(ctx, p);
  }

  template <>
  std::optional<std::vector<int>> convertToVector<int>(
      Context& ctx, const Parameter& p) noexcept {
    return convertToVector_impl<int>(ctx, p);
  }

}  // namespace mfem_mgis::internals

namespace mfem_mgis {

  InvalidResult Parameter::reportUnmatchedParameterType(Context& ctx) noexcept {
    return ctx.registerErrorMessage(
        "the type of parameter is not the expected one");
  }  // end of reportUnmatchedParameterType

  void Parameter::raiseUnmatchedParameterType(attributes::Throwing) {
    raise(
        "Parameter::raiseUnmatchedParameterType: "
        "the type of parameter is not the expected one");
  }  // end of raiseUnmatchedParameterType

  Parameter::Parameter() = default;
  Parameter::Parameter(Parameter&&) = default;
  Parameter::Parameter(const Parameter&) = default;
  Parameter& Parameter::operator=(Parameter&&) = default;
  Parameter& Parameter::operator=(const Parameter&) = default;

  Parameter::Parameter(const char* const src)
      : ParameterVariant(std::string(src)) {}  // end of Parameter

  Parameter::Parameter(std::string_view src)
      : ParameterVariant(std::string(src)) {}  // end of Parameter

  Parameter& Parameter::operator=(const char* const src) {
    this->operator=(std::string(src));
    return *this;
  }  // end of Parameter::operator=

  Parameter& Parameter::operator=(std::string_view src) {
    this->operator=(std::string(src));
    return *this;
  }  // end of Parameter::operator=

  Parameter::~Parameter() = default;

  Parameter get(attributes::Throwing, const Parameters& p, std::string_view n) {
    return p.get(throwing, n);
  }  // end of get

  OptionalReference<const Parameter> get(Context& ctx,
                                         const Parameters& p,
                                         std::string_view n) noexcept {
    return p.get(ctx, n);
  }  // end of get

  Parameter get_if(const Parameters& p,
                   std::string_view n,
                   const Parameter& v) noexcept {
    if (contains(p, n)) {
      return p.get(throwing, n);
    }
    return v;
  }  // end of get

}  // end of namespace mfem_mgis
