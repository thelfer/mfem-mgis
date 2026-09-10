/*!
 * \file   MFEMMGIS/Parameter.ixx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETER_IXX
#define LIB_MFEM_MGIS_PARAMETER_IXX

namespace mfem_mgis {

  inline ParameterVariant& Parameter::as_std_variant() noexcept {
    return *this;
  }  // end of as_std_variant

  inline const ParameterVariant& Parameter::as_std_variant() const noexcept {
    return *this;
  }  // end of as_std_variant

  template <typename ResultType>
  bool is(const Parameter& p) noexcept {
    return std::holds_alternative<ResultType>(p.as_std_variant());
  }  // end of is

  template <>
  inline bool is<double>(const Parameter& p) noexcept {
    return std::holds_alternative<double>(p.as_std_variant()) ||
           std::holds_alternative<int>(p.as_std_variant());
  }

  template <typename ResultType>
  OptionalGetResultType<ResultType> get(Context& ctx,
                                        const Parameter& p) noexcept {
    if (!is<ResultType>(p)) {
      return Parameter::reportUnmatchedParameterType(ctx);
    }
    return OptionalReference<const ResultType>(
        &(std::get<ResultType>(p.as_std_variant())));
  }  // end of get

  template <>
  inline OptionalGetResultType<double> get<double>(
      Context& ctx, const Parameter& p) noexcept {
    if (!is<double>(p)) {
      return Parameter::reportUnmatchedParameterType(ctx);
    }
    if (std::holds_alternative<int>(p.as_std_variant())) {
      return std::get<int>(p.as_std_variant());
    }
    return std::get<double>(p.as_std_variant());
  }  // end of get

  template <typename ResultType>
  GetResultType<ResultType> get(attributes::Throwing, const Parameter& p) {
    if (!is<ResultType>(p)) {
      Parameter::raiseUnmatchedParameterType(throwing);
    }
    return std::get<ResultType>(p.as_std_variant());
  }  // end of get

  template <>
  inline GetResultType<double> get<double>(attributes::Throwing,
                                           const Parameter& p) {
    if (!is<double>(p)) {
      Parameter::raiseUnmatchedParameterType(throwing);
    }
    if (std::holds_alternative<int>(p.as_std_variant())) {
      return std::get<int>(p.as_std_variant());
    }
    return std::get<double>(p.as_std_variant());
  }  // end of get

  template <typename ResultType>
  bool is(attributes::Throwing, const Parameters& p, std::string_view n) {
    const auto& v = p.get(throwing, n);
    return is<ResultType>(v);
  }  // end of is

  template <typename ResultType>
  OptionalGetResultType<ResultType> get(Context& ctx,
                                        const Parameters& p,
                                        std::string_view n) noexcept {
    const auto& ov = p.get(ctx, n);
    if (isInvalid(ov)) {
      return {};
    }
    if (!is<ResultType>(*ov)) {
      return Parameters::reportUnmatchedParameterType(ctx, n);
    }
    return get<ResultType>(ctx, *ov);
  }  // end of get

  template <>
  inline OptionalGetResultType<double> get<double>(
      Context& ctx, const Parameters& p, std::string_view n) noexcept {
    const auto& ov = p.get(ctx, n);
    if (isInvalid(ov)) {
      return {};
    }
    if (!is<double>(*ov)) {
      return Parameters::reportUnmatchedParameterType(ctx, n);
    }
    return get<double>(ctx, *ov);
  }  // end of get

  template <typename ResultType>
  std::optional<ResultType> get_if(Context& ctx,
                                   const Parameters& p,
                                   std::string_view n,
                                   const ResultType& v) noexcept {
    if (contains(p, n)) {
      const auto ovalue = get<ResultType>(ctx, p, n);
      if (isInvalid(ovalue)) {
        return {};
      }
      return *ovalue;
    }
    return v;
  }  // end of get

  template <typename ResultType>
  GetResultType<ResultType> get(attributes::Throwing,
                                const Parameters& p,
                                std::string_view n) {
    const auto& v = p.get(throwing, n);
    if (!is<ResultType>(v)) {
      Parameters::raiseUnmatchedParameterType(throwing, n);
    }
    return get<ResultType>(throwing, v);
  }  // end of get

  template <>
  inline GetResultType<double> get<double>(attributes::Throwing,
                                           const Parameters& p,
                                           std::string_view n) {
    const auto& v = p.get(throwing, n);
    if (!is<double>(v)) {
      Parameters::raiseUnmatchedParameterType(throwing, n);
    }
    return get<double>(throwing, v);
  }  // end of get

  template <typename ResultType>
  ResultType get_if(attributes::Throwing,
                    const Parameters& p,
                    std::string_view n,
                    const ResultType& v) {
    if (contains(p, n)) {
      return get<ResultType>(throwing, p, n);
    }
    return v;
  }  // end of get

  namespace internals {

    template <typename ValueType>
    requires((std::same_as<ValueType, real>) || (std::same_as<ValueType, int>))
        [[nodiscard]] std::optional<std::vector<ValueType>> convertToVector(
            Context&, const Parameter&) noexcept;

    template <>
    MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<real>>
    convertToVector<real>(Context&, const Parameter&) noexcept;

    template <>
    MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<int>>
    convertToVector<int>(Context&, const Parameter&) noexcept;

  }  // end of namespace internals

  template <typename ValueType>
  [[nodiscard]] auto convert(Context& ctx, const Parameter& p) noexcept {
    if constexpr (std::same_as<ValueType, std::vector<real>>) {
      return ::mfem_mgis::internals::convertToVector<real>(ctx, p);
    } else if constexpr (std::same_as<ValueType, std::vector<int>>) {
      return ::mfem_mgis::internals::convertToVector<int>(ctx, p);
    } else {
      return get<ValueType>(ctx, p);
    }
  }  // end of convert

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARAMETER_IXX */
