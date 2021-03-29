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
  static Parameters& insert_implementation(Parameters& parameters,
                                           const Container& src) {
    for (const auto& p : src) {
      parameters.insert(p.first, p.second);
    }
    return parameters;
  }  // end of insert_implementation

  void Parameters::raiseUnmatchedParameterType(std::string_view n) {
    std::string msg("Parameters::raiseUnmatchedParameterType :");
    msg += "the type of parameter '";
    msg += n;
    msg += "' is not the expected one";
    mgis::raise(msg);
  }  // end of raiseUnmatchedParameterType

  Parameters::Parameters() = default;

  Parameters::Parameters(const Parameters&) = default;

  Parameters::Parameters(Parameters&&) = default;

  Parameters& Parameters::operator=(const Parameters&) = default;

  Parameters& Parameters::operator=(Parameters&&) = default;

  Parameters::const_iterator Parameters::begin() const {
    return std::map<std::string, Parameter, std::less<>>::cbegin();
  }  // end of begin

  Parameters::const_iterator Parameters::cbegin() const {
    return std::map<std::string, Parameter, std::less<>>::cbegin();
  }  // end of cbegin

  Parameters::const_iterator Parameters::end() const {
    return std::map<std::string, Parameter, std::less<>>::cend();
  }  // end of end

  Parameters::const_iterator Parameters::cend() const {
    return std::map<std::string, Parameter, std::less<>>::cend();
  }  // end of cend

  bool Parameters::contains(std::string_view n) const {
    return this->count(n) != 0;
  }  // end of contains

  Parameters& Parameters::insert(const Parameters& src) {
    return insert_implementation(*this, src);
  }  // end of insert

  Parameters& Parameters::insert(const std::map<std::string, Parameter>& src) {
    return insert_implementation(*this, src);
  }  // end of insert

  Parameters& Parameters::insert(
      const std::initializer_list<std::map<std::string, Parameter>::value_type>&
          src) {
    return insert_implementation(*this, src);
  }  // end of insert

  Parameters& Parameters::insert(std::string_view n, const Parameter& p) {
    if (this->count(n) != 0) {
      std::string msg("Parameters::insert: parameter '");
      msg += n;
      msg += "' has already been declared";
    }
    std::map<std::string, Parameter, std::less<>>::value_type v{n, p};
    std::map<std::string, Parameter, std::less<>>::insert(std::move(v));
    return *this;
  }  // end of insert

  const Parameter& Parameters::get(std::string_view n) const {
    const auto i = this->find(n);
    if (i == this->end()) {
      std::string msg("Parameters::get: parameter '");
      msg += n;
      msg += "' is not declared";
      mgis::raise(msg);
    }
    return i->second;
  }

  Parameters::~Parameters() = default;

  bool contains(const Parameters& p, std::string_view n) {
    return p.contains(n);
  }  // end of contains

}  // end of namespace mfem_mgis
