/*!
 * \file   MFEMMGIS/Parameter.ixx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETER_IXX
#define LIB_MFEM_MGIS_PARAMETER_IXX

namespace mfem_mgis {

  template <typename ResultType>
  bool is(const Parameter& p) {
    return std::holds_alternative<ResultType>(p);
  }  // end of is

  template <>
  inline bool is<double>(const Parameter& p) {
    return std::holds_alternative<double>(p) || std::holds_alternative<int>(p);
  }

  template <typename ResultType>
  GetResultType<ResultType> get(const Parameter& p) {
    if (!is<ResultType>(p)) {
      Parameter::raiseUnmatchedParameterType();
    }
    return std::get<ResultType>(p);
  }  // end of get

  template <>
  inline GetResultType<double> get<double>(const Parameter& p) {
    if ((!is<double>(p)) && (!is<int>(p))) {
      Parameter::raiseUnmatchedParameterType();
    }
    if (std::holds_alternative<int>(p)) {
      return std::get<int>(p);
    }
    return std::get<double>(p);
  }  // end of get

  template <typename ResultType>
  bool is(const Parameters& p, std::string_view n) {
    const auto& v = p.get(n);
    return is<ResultType>(v);
  }  // end of is

  template <typename ResultType>
  GetResultType<ResultType> get(const Parameters& p, std::string_view n) {
    const auto& v = p.get(n);
    if (!is<ResultType>(v)) {
      Parameters::raiseUnmatchedParameterType(n);
    }
    return get<ResultType>(v);
  }  // end of get

  template <>
  inline GetResultType<double> get<double>(const Parameters& p,
                                           std::string_view n) {
    const auto& v = p.get(n);
    if (!is<double>(v)) {
      Parameters::raiseUnmatchedParameterType(n);
    }
    return get<double>(v);
  }  // end of get

  template <typename ResultType>
  ResultType get_if(const Parameters& p,
                    std::string_view n,
                    const ResultType& v) {
    if (contains(p, n)) {
      return get<ResultType>(p, n);
    }
    return v;
  }  // end of get

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARAMETER_IXX */
