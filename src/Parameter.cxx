/*!
 * \file   src/Parameter.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/03/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameter.hxx"

namespace mfem_mgis {

  void Parameter::raiseUnmatchedParameterType() {
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

  Parameter get(const Parameters& p, std::string_view n) {
    return p.get(n);
  }  // end of get

  Parameter get_if(const Parameters& p,
                   std::string_view n,
                   const Parameter& v) {
    if (contains(p, n)) {
      return p.get(n);
    }
    return v;
  }  // end of get

}  // end of namespace mfem_mgis