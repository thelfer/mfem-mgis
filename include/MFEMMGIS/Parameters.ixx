/*!
 * \file   MFEMMGIS/Parameters.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETERS_IXX
#define LIB_MFEM_MGIS_PARAMETERS_IXX

namespace mfem_mgis {

  template <typename ResultType>
  bool is(const Parameters& p, std::string_view n) {
    const auto& v = p.get(n);
    return std::holds_alternative<ResultType>(v);
  }  // end of is

  template <typename ResultType>
  const ResultType& get(const Parameters& p, std::string_view n) {
    const auto& v = p.get(n);
    if (!std::holds_alternative<ResultType>(v)) {
      Parameters::raiseUnmatchedParameterType(n);
    }
    return std::get<ResultType>(v);
  }  // end of get

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARAMETERS_IXX */
