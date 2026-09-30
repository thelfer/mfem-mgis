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

  //! \brief concept satisfied by the types of the values of a parameter
  template <typename T>
  concept ParameterValueConcept =
      ((std::same_as<T, bool>) || (std::same_as<T, size_type>) ||
       (std::same_as<T, real>) || (std::same_as<T, std::string>) ||
       (std::same_as<T, std::vector<Parameter>>) ||
       (std::same_as<T, Parameters>) ||
       (std::same_as<T, std::function<real(const real)>>));

  /* aliases to equivalent python' types */

  //! \brief equivalent of the python list
  using list = std::vector<Parameter>;
  //! \brief equivalent of the python dict
  using dict = Parameters;

  /*!
   * \brief variant class used to initialize objects.
   * \note this class wraps a `ParameterVariant` to provide a convenient way to
   * handle different types of parameters.
   */
  struct MFEM_MGIS_EXPORT [[nodiscard]] Parameter : private ParameterVariant {
    /*!
     * \brief report that the type of a parameter is not the expected one.
     * \param[in, out] ctx: execution context
     * \return an invalid result
     */
    static InvalidResult reportUnmatchedParameterType(Context& ctx) noexcept;
    /*!
     * \brief throw an exception stating that the parameter type is not the
     * expected one.
     * \param[in] throwing: dummy attribute to indicate that this function may
     * throw an exception
     */
    [[noreturn]] static void raiseUnmatchedParameterType(
        attributes::Throwing throwing);
    /*!
     * \brief create a parameter holding a `std::vector<Parameter>` containing a
     * copy of the given vector.
     * \param[in] values: values to be inserted
     * \return the created parameter
     */
    template <ParameterValueConcept ParameterType>
    [[nodiscard]] static Parameter from(
        const std::vector<ParameterType>& values) noexcept;
    /*!
     * \brief create a parameter holding a `std::vector<Parameter>` containing a
     * copy of the given set.
     * \param[in] values: values to be inserted
     * \return the created parameter
     */
    template <ParameterValueConcept ParameterType>
    [[nodiscard]] static Parameter from(
        const std::set<ParameterType>& values) noexcept;
    // inheriting constructors
    using ParameterVariant::ParameterVariant;
    //! \brief default constructor
    Parameter();
    //! \brief move constructor
    Parameter(Parameter&&) noexcept;
    //! \brief copy constructor
    Parameter(const Parameter&);
    /*!
     * \brief constructor from a C-string
     * \param[in] src: source
     */
    Parameter(const char* const src);
    /*!
     * \brief constructor from a `std::string_view`
     * \param[in] src: source
     */
    Parameter(std::string_view src);
    //! \brief move assignment
    Parameter& operator=(Parameter&&) noexcept;
    //! \brief copy assignment
    Parameter& operator=(const Parameter&);
    //! \brief inheriting assignment operators
    using ParameterVariant::operator=;
    /*!
     * \brief assignment from a C-string
     * \param[in] src: source
     * \return a reference to this object
     */
    Parameter& operator=(const char* const src);
    /*!
     * \brief assignment from a `std::string_view`
     * \param[in] src: source
     * \return a reference to this object
     */
    Parameter& operator=(std::string_view src);
    //! \return the underlying variant
    ParameterVariant& as_std_variant() noexcept;
    //! \return the underlying variant
    const ParameterVariant& as_std_variant() const noexcept;
    //! \brief destructor
    ~Parameter();
  };  // end of struct Parameter

  //! \brief Type alias for result type
  template <typename ResultType>
  using GetResultType = std::conditional_t<std::is_same_v<ResultType, double>,
                                           ResultType,
                                           const ResultType&>;
  //! \brief Type alias for optional result type
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
  [[nodiscard]] bool is(const Parameter& p) noexcept;

  /*!
   * \brief specialisation of the `is` function for double
   * \param[in] p: parameter
   * \return true if the given parameter is a double or an integer
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
      Context& ctx, const Parameter& p) noexcept;

  /*!
   * \brief specialisation of the `get` function for double
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter
   * \return the value of the parameter if it is a double or an integer
   */
  template <>
  [[nodiscard]] OptionalGetResultType<double> get<double>(
      Context&, const Parameter&) noexcept;
  /*!
   * \brief get the value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameter
   * \return the value of the parameter
   * \throws std::runtime_error if the parameter does not have the expected
   * type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing throwing,
                                              const Parameter& p);

  /*!
   * \brief specialisation of the `get` function for double
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameter
   * \return the value of the parameter
   * \throws std::runtime_error if the parameter is neither a double nor an
   * integer.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing throwing,
                                                  const Parameter&);

  /*!
   * \brief check if the given parameter exists
   * \param[in] p: parameters
   * \param[in] n: name
   * \return true if the given parameter exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool contains(const Parameters& p,
                                               std::string_view n) noexcept;

  /*!
   * \brief check if the given parameter has the given type
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return true if the given parameter has the given type
   * \throws std::runtime_error if the parameter does not exist
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] bool is(attributes::Throwing throwing,
                        const Parameters& p,
                        std::string_view n);
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
      Context& ctx, const Parameters& p, std::string_view n) noexcept;
  /*!
   * \brief get the value of the parameter if it exists
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter if it exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] OptionalReference<const Parameter> get(
      Context& ctx, const Parameters& p, std::string_view n) noexcept;
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
   * \return value of the parameter if present, a default value otherwise.
   * Empty on failure.
   */
  template <typename ResultType>
  [[nodiscard]] std::optional<ResultType> get_if(Context& ctx,
                                                 const Parameters& p,
                                                 std::string_view n,
                                                 const ResultType& v) noexcept;
  /*!
   * \brief get the value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws std::runtime_error if the parameter does not exist or does not have
   * the expected type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing throwing,
                                              const Parameters& p,
                                              std::string_view n);
  /*!
   * \brief get the value of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws std::runtime_error if the parameter does not exist
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get(attributes::Throwing throwing,
                                               const Parameters& p,
                                               std::string_view n);

  /*!
   * \brief get the value of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \return the value of the parameter
   * \throws std::runtime_error if the parameter does not exist or is not a
   * number
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing throwing,
                                                  const Parameters&,
                                                  std::string_view);

  /*!
   * \brief get the value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \return the value of the parameter if present, a default value otherwise
   * \throws std::runtime_error if the parameter exists but does not have the
   * expected type.
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  template <typename ResultType>
  [[nodiscard]] ResultType get_if(attributes::Throwing throwing,
                                  const Parameters& p,
                                  std::string_view n,
                                  const ResultType& v);
  /*!
   * \brief get the value of the parameter if present, a default value otherwise
   * \return the value of the parameter if present, a default value otherwise
   *
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get_if(const Parameters& p,
                                                  std::string_view n,
                                                  const Parameter& v) noexcept;
  /*!
   * \brief convert the parameter to the given type
   * \tparam ValueType: target type
   * \param[in, out] ctx: execution context
   * \param[in] p: parameter
   * \return the parameter converted to the given type, empty on failure
   */
  template <typename ValueType>
  [[nodiscard]] auto convert(Context& ctx, const Parameter& p) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.ixx"
#include "MFEMMGIS/Utilities/ParametersValidator.hxx"

#endif /* LIB_MFEM_MGIS_PARAMETER_HXX */
