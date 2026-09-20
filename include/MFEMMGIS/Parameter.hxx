/*!
 * \file   include/MFEMMGIS/Parameter.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETER_HXX
#define LIB_MFEM_MGIS_PARAMETER_HXX

#include <set>
#include <variant>
#include <concepts>
#include <functional>
#include <type_traits>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;

  //! \brief variant type for parameter values
  using ParameterVariant = std::variant<std::monostate,
                                        bool,
                                        size_type,
                                        real,
                                        std::string,
                                        std::vector<Parameter>,
                                        Parameters,
                                        std::function<real(const real)>>;

  template <typename T>
  concept ParameterValueConcept =
      ((std::same_as<T, bool>) || (std::same_as<T, size_type>) ||
       (std::same_as<T, real>) || (std::same_as<T, std::string>) ||
       (std::same_as<T, std::vector<Parameter>>) ||
       (std::same_as<T, Parameters>) ||
       (std::same_as<T, std::function<real(const real)>>));

  /* aliases to equivalent python' types */

  using list = std::vector<Parameter>;
  using dict = Parameters;

  /*!
   * \brief variant class used to initialize objects.
   * \note this class wraps a `ParameterVariant` to provide a convenient way to
   * handle different types of parameters.
   */
  struct MFEM_MGIS_EXPORT [[nodiscard]] Parameter : private ParameterVariant {
    /*!
     * \brief report that the type of a parameter is the expected one.
     * \param[in, out] ctx: execution context
     */
    static InvalidResult reportUnmatchedParameterType(Context&) noexcept;
    /*!
     * \brief throw an exception if the parameter type is not the expected one.
     */
    [[noreturn]] static void raiseUnmatchedParameterType(attributes::Throwing);
    /*!
     * \brief create a parameter holding a `std::vector<Parameter>` containing a
     * copy of the given vector.
     * \param[in] values: values to be inserted
     */
    template <ParameterValueConcept ParameterType>
    [[nodiscard]] static Parameter from(
        const std::vector<ParameterType>&) noexcept;
    /*!
     * \brief create a parameter holding a `std::vector<Parameter>` containing a
     * copy of the given set.
     * \param[in] values: values to be inserted
     */
    template <ParameterValueConcept ParameterType>
    [[nodiscard]] static Parameter from(
        const std::set<ParameterType>&) noexcept;
    // inheriting constructors
    using ParameterVariant::ParameterVariant;
    //! \brief default constructor
    Parameter();
    //! \brief move constructor
    Parameter(Parameter&&);
    //! \brief copy constructor
    Parameter(const Parameter&);
    /*!
     * \brief constructor from a C-string
     * \param[in] src: source
     */
    Parameter(const char* const);
    /*!
     * \brief constructor from a `std::string_view`
     * \param[in] src: source
     */
    Parameter(std::string_view);
    //! \brief move assignment
    Parameter& operator=(Parameter&&);
    //! \brief copy assignment
    Parameter& operator=(const Parameter&);
    //! \brief inheriting assignment operators
    using ParameterVariant::operator=;
    /*!
     * \brief assignment from a C-string
     * \param[in] src: source
     */
    Parameter& operator=(const char* const);
    /*!
     * \brief assignment from a `std::string_view`
     * \param[in] src: source
     */
    Parameter& operator=(std::string_view);
    //
    ParameterVariant& as_std_variant() noexcept;
    const ParameterVariant& as_std_variant() const noexcept;
    //! \brief destructor
    ~Parameter();
  };  // end of struct Parameter

  //! \brief Type alias for result type
  template <typename ResultType>
  using GetResultType = std::conditional_t<std::is_same_v<ResultType, double>,
                                           ResultType,
                                           const ResultType&>;
  template <typename ResultType>
  using OptionalGetResultType =
      std::conditional_t<std::is_same_v<ResultType, double>,
                         std::optional<ResultType>,
                         OptionalReference<const ResultType>>;

  /*!
   * \brief check if the given parameter has the given type
   * \param[in] p: parameter
   * \return true if the given parameter has the given type
   */
  template <typename ResultType>
  [[nodiscard]] bool is(const Parameter&) noexcept;

  /*!
   * \brief partial specialisation of the `is` function for double
   * \param[in] p: parameter
   * \return true if the given parameter is a double
   */
  template <>
  [[nodiscard]] bool is<double>(const Parameter&) noexcept;
  /*!
   * \brief get the value of the parameter if it has the expected type
   * \tparam ResultType: expected type of the parameter
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter
   * \return the value of the parameter if it has the expected type
   */
  template <typename ResultType>
  [[nodiscard]] OptionalGetResultType<ResultType> get(
      Context&, const Parameter&) noexcept;

  /*!
   * \brief partial specialisation of the `get` function for double
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter
   * \return the value of the parameter if it is a double
   */
  template <>
  [[nodiscard]] OptionalGetResultType<double> get<double>(
      Context&, const Parameter&) noexcept;
  /*!
   * \brief get the value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameter
   * \return the value of the parameter
   * \throws if the parameter does not have the expected type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing,
                                              const Parameter&);

  /*!
   * \brief partial specialisation of the `get` function for double
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameter
   * \return the value of the parameter
   * \throws if the parameter is not a double.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing,
                                                  const Parameter&);

  /*!
   * \brief check if the given parameter exists
   * \param[in] p: parameters
   * \param[in] n: name
   * \return true if the given parameter exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool contains(const Parameters&,
                                               std::string_view) noexcept;

  /*!
   * \brief check if the given parameter has the given type
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return true if the given parameter has the given type
   * \throws if the parameter does not exist
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] bool is(attributes::Throwing,
                        const Parameters&,
                        std::string_view);
  /*!
   * \brief get the value of the parameter if it exists and has the expected
   * type
   * \tparam ResultType: expected type of the parameter
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter if it exists and has the expected type
   */
  template <typename ResultType>
  [[nodiscard]] OptionalGetResultType<ResultType> get(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \brief get the value of the parameter if it exists
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter if it exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] OptionalReference<const Parameter> get(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \brief get the value of the parameter if it exists and is a number
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter if it exists and is a number
   */
  template <>
  [[nodiscard]] OptionalGetResultType<double> get<double>(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \brief get the value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \return value of the parameter if present, a default value otherwise
   */
  template <typename ResultType>
  [[nodiscard]] std::optional<ResultType> get_if(Context&,
                                                 const Parameters&,
                                                 std::string_view,
                                                 const ResultType&) noexcept;
  /*!
   * \brief get the value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws if the parameter does not exist or does not have the expected type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing,
                                              const Parameters&,
                                              std::string_view);
  /*!
   * \brief get the value of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws if the parameter does not exist
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get(attributes::Throwing,
                                               const Parameters&,
                                               std::string_view);

  /*!
   * \brief get the value of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws if the parameter does not exist or is not a number
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing,
                                                  const Parameters&,
                                                  std::string_view);

  /*!
   * \brief get the value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \return the value of the parameter if present, a default value otherwise
   * \throws if the parameter exists but does not have the expected type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] ResultType get_if(attributes::Throwing,
                                  const Parameters&,
                                  std::string_view,
                                  const ResultType&);
  /*!
   * \return the value of the parameter if present, a default value otherwise
   *
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get_if(const Parameters&,
                                                  std::string_view,
                                                  const Parameter&) noexcept;
  /*!
   * \brief convert the parameter to the given type
   * \tparam ValueType: target type
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter
   * \return he parameter converted to the given type
   */
  template <typename ValueType>
  [[nodiscard]] auto convert(Context&, const Parameter&) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.ixx"
#include "MFEMMGIS/Utilities/ParametersValidator.hxx"

#endif /* LIB_MFEM_MGIS_PARAMETER_HXX */
