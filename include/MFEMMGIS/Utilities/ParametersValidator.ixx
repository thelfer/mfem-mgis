/*!
 * \file   include/MFEMMGIS/Utilities/ParametersValidator.ixx
 * \brief  This files implements the inline methods of the `ParametersValidator`
 * class
 * \author Thomas Helfer
 * \date   19/09/2026
 */

#ifndef LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_IXX
#define LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_IXX

namespace mfem_mgis {

  template <typename... Types>
  ParametersValidator& ParametersValidator::add(
      const std::string& k, const AddArguments& opts) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>))) {
    auto v = [k](Context& ctx, const Parameter& d) noexcept -> bool {
      const auto b = (... || is<Types>(d));
      if (!b) {
        if constexpr (sizeof...(Types) == 1) {
          return ParametersValidator::reportUnmatchedTypeError(ctx, k);
        } else {
          return ParametersValidator::reportUnmatchedTypesError(ctx, k);
        }
      }
      return true;
    };
    this->add(k, v, opts);
    return *this;
  }  // end of add

  template <typename... Types>
  ParametersValidator& ParametersValidator::add(
      const std::string& k,
      const std::string& d,
      const AddArguments& opts) noexcept
      requires((sizeof...(Types) > 0) &&
               (... && (ParameterValueConcept<std::decay_t<Types>>))) {
    this->add(k, d, opts);  // pre-declare the key with the documentation
    return this->template add<Types...>(k, opts);
  }  // end of add

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_UTILITIES_PARAMETERSVALIDATOR_IXX */
