/*!
 * \file   src/Parameter.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  template <typename Container>
  static Parameters& insert_implementation(attributes::Throwing,
                                           Parameters& parameters,
                                           const Container& src) {
    for (const auto& p : src) {
      parameters.insert(throwing, p.first, p.second);
    }
    return parameters;
  }  // end of insert_implementation

  InvalidResult Parameters::reportMissingKey(Context& ctx,
                                             std::string_view n) noexcept {
    auto msg = std::string{"no parameter '"};
    msg += n;
    msg += "' declared";
    return ctx.registerErrorMessage(msg);
  }  // end of reportMissingKey

  InvalidResult Parameters::reportUnmatchedParameterType(
      Context& ctx, std::string_view n) noexcept {
    auto msg = std::string{"the type of parameter '"};
    msg += n;
    msg += "' is not the expected one";
    return ctx.registerErrorMessage(msg);
  }

  void Parameters::raiseUnmatchedParameterType(attributes::Throwing,
                                               std::string_view n) {
    std::string msg("Parameters::raiseUnmatchedParameterType :");
    msg += "the type of parameter '";
    msg += n;
    msg += "' is not the expected one";
    raise(msg);
  }  // end of raiseUnmatchedParameterType

  Parameters::Parameters() noexcept = default;

  Parameters::Parameters(const Parameters&) noexcept = default;

  Parameters::Parameters(Parameters&&) noexcept = default;

  Parameters& Parameters::operator=(const Parameters&) noexcept = default;

  Parameters& Parameters::operator=(Parameters&&) noexcept = default;

  Parameters::const_iterator Parameters::begin() const noexcept {
    return std::map<std::string, Parameter, std::less<>>::cbegin();
  }  // end of begin

  Parameters::const_iterator Parameters::cbegin() const noexcept {
    return std::map<std::string, Parameter, std::less<>>::cbegin();
  }  // end of cbegin

  Parameters::const_iterator Parameters::end() const noexcept {
    return std::map<std::string, Parameter, std::less<>>::cend();
  }  // end of end

  Parameters::const_iterator Parameters::cend() const noexcept {
    return std::map<std::string, Parameter, std::less<>>::cend();
  }  // end of cend

  bool Parameters::contains(std::string_view n) const noexcept {
    return this->count(n) != 0;
  }  // end of contains

  Parameters& Parameters::insert(attributes::Throwing, const Parameters& src) {
    return insert_implementation(throwing, *this, src);
  }  // end of insert

  Parameters& Parameters::insert(attributes::Throwing,
                                 const std::map<std::string, Parameter>& src) {
    return insert_implementation(throwing, *this, src);
  }  // end of insert

  Parameters& Parameters::insert(
      attributes::Throwing,
      const std::initializer_list<std::map<std::string, Parameter>::value_type>&
          src) {
    return insert_implementation(throwing, *this, src);
  }  // end of insert

  Parameters& Parameters::insert(attributes::Throwing,
                                 std::string_view n,
                                 const Parameter& p) {
    if (this->count(n) != 0) {
      std::string msg("Parameters::insert: parameter '");
      msg += n;
      msg += "' has already been declared";
      raise(msg);
    }
    std::map<std::string, Parameter, std::less<>>::value_type v{n, p};
    std::map<std::string, Parameter, std::less<>>::insert(std::move(v));
    return *this;
  }  // end of insert

  bool Parameters::insert(Context& ctx,
                          std::string_view n,
                          const Parameter& p) noexcept {
    if (this->count(n) != 0) {
      return ctx.registerErrorMessage("parameter '" + std::string{n} +
                                      "' has already been declared");
    }
    std::map<std::string, Parameter, std::less<>>::value_type v{n, p};
    std::map<std::string, Parameter, std::less<>>::insert(std::move(v));
    return true;
  }  // end of insert

  OptionalReference<const Parameter> Parameters::get(
      Context& ctx, std::string_view n) const noexcept {
    const auto i = this->find(n);
    if (i == this->end()) {
      return ctx.registerErrorMessage("parameter '" + std::string{n} +
                                      "' is not declared");
    }
    return OptionalReference<const Parameter>(&(i->second));
  }  // end of get

  const Parameter& Parameters::get(attributes::Throwing,
                                   std::string_view n) const {
    const auto i = this->find(n);
    if (i == this->end()) {
      std::string msg("Parameters::get: parameter '");
      msg += n;
      msg += "' is not declared";
      raise(msg);
    }
    return i->second;
  }

  Parameters::~Parameters() = default;

  bool contains(const Parameters& p, std::string_view n) noexcept {
    return p.contains(n);
  }  // end of contains

  bool checkParameters(Context& ctx,
                       const Parameters& parameters,
                       const std::vector<std::string>& names) noexcept {
    for (const auto& p : parameters) {
      if (std::find(names.begin(), names.end(), p.first) == names.end()) {
        auto msg = std::string("checkParameters: invalid parameter '" +
                               p.first + "'.");
        if (!names.empty()) {
          msg += "\nAllowed parameters are:\n";
          for (const auto& n : names) {
            msg += "- " + n + '\n';
          }
        } else {
          msg += "No parameters allowed";
        }
        return ctx.registerErrorMessage(msg);
      }
    }
    return true;
  }  // end of checkParameters

  bool checkParameters(
      Context& ctx,
      const Parameters& parameters,
      const std::map<std::string, std::string>& descriptions) noexcept {
    for (const auto& p : parameters) {
      if (!descriptions.contains(p.first)) {
        auto msg = std::string("checkParameters: invalid parameter '" +
                               p.first + "'.");
        if (!descriptions.empty()) {
          msg += "\nAllowed parameters are:\n";
          for (const auto& [n, d] : descriptions) {
            msg += "- " + n + ": " + d + '\n';
          }
        } else {
          msg += "No parameters allowed";
        }
        return ctx.registerErrorMessage(msg);
      }
    }
    return true;
  }  // end of checkParameters

  void checkParameters(attributes::Throwing,
                       const Parameters& parameters,
                       const std::vector<std::string>& names) {
    for (const auto& p : parameters) {
      if (std::find(names.begin(), names.end(), p.first) == names.end()) {
        auto msg = std::string("checkParameters: invalid parameter '" +
                               p.first + "'.");
        if (!names.empty()) {
          msg += "\nAllowed parameters are:\n";
          for (const auto& n : names) {
            msg += "- " + n + '\n';
          }
        } else {
          msg += "No parameters allowed";
        }
        raise(msg);
      }
    }
  }  // end of checkParameters

  void checkParameters(attributes::Throwing,
                       const Parameters& parameters,
                       const std::map<std::string, std::string>& descriptions) {
    for (const auto& p : parameters) {
      if (!descriptions.contains(p.first)) {
        auto msg = std::string("checkParameters: invalid parameter '" +
                               p.first + "'.");
        if (!descriptions.empty()) {
          msg += "\nAllowed parameters are:\n";
          for (const auto& [n, d] : descriptions) {
            msg += "- " + n + ": " + d + '\n';
          }
        } else {
          msg += "No parameters allowed";
        }
        raise(msg);
      }
    }
  }  // end of checkParameters

  Parameters extract(attributes::Throwing,
                     const Parameters& parameters,
                     const std::vector<std::string>& names) {
    auto r = Parameters{};
    for (const auto& n : names) {
      if (contains(parameters, n)) {
        r.insert(throwing, n, parameters.get(throwing, n));
      }
    }
    return r;
  }  // end of extract

  std::optional<Parameters> extract(
      Context& ctx,
      const Parameters& parameters,
      const std::vector<std::string>& names) noexcept {
    auto r = Parameters{};
    for (const auto& n : names) {
      if (contains(parameters, n)) {
        const auto ovalue = parameters.get(ctx, n);
        if (isInvalid(ovalue)) {
          return {};
        }
        if (!r.insert(ctx, n, *ovalue)) {
          return {};
        }
      }
    }
    return r;
  }  // end of extract

  [[nodiscard]] static std::optional<std::pair<std::string, Parameters>>
  extractFactoryArgumentImplementation(Context& ctx,
                                       const Parameters& parameters) noexcept {
    if (parameters.size() != 1u) {
      return ctx.registerErrorMessage(
          "the parameter must be a dictionary with a unique entry");
    }
    const auto& [n, params] = *(parameters.begin());
    if (!is<Parameters>(params)) {
      return ctx.registerErrorMessage(
          "expected a dictionary to define the parameter of the '" + n +
          "' object");
    }
    return std::pair<std::string, Parameters>{
        n, get<Parameters>(throwing, params)};
  }  // end of extractFactoryArgumentImplementation

  std::optional<std::pair<std::string, Parameters>> extractFactoryArgument(
      Context& ctx, const Parameters& parameters) noexcept {
    return extractFactoryArgumentImplementation(ctx, parameters);
  }  // end of extractFactoryArgument

  std::pair<std::string, Parameters> extractFactoryArgument(
      attributes::Throwing, const Parameters& parameters) {
    auto ctx = Context{};
    const auto ovalue = extractFactoryArgumentImplementation(ctx, parameters);
    if (isInvalid(ovalue)) {
      raise(ctx.getErrorMessage());
    }
    return *ovalue;
  }  // end of extractFactoryArgument

}  // end of namespace mfem_mgis
