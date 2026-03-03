/*!
 * \file   include/MFEMMGIS/Parameter.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETER_HXX
#define LIB_MFEM_MGIS_PARAMETER_HXX

#include <variant>
#include <functional>
#include <type_traits>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;

  //! \brief a simple alias
  using ParameterVariant = std::variant<std::monostate,
                                        bool,
                                        size_type,
                                        real,
                                        std::string,
                                        std::vector<Parameter>,
                                        std::vector<int>,
                                        Parameters,
                                        std::function<real(const real)>>;

  /*!
   *
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
    // inheriting constructors
    using ParameterVariant::ParameterVariant;
    // \brief default constructor
    Parameter();
    // \brief move constructor
    Parameter(Parameter&&);
    // \brief copy constructor
    Parameter(const Parameter&);
    /*!
     * \brief constructor from C-string
     * \param[in] src: source
     */
    Parameter(const char* const);
    /*!
     * \brief constructor from a string_view
     * \param[in] src: source
     */
    Parameter(std::string_view);
    // \brief move assignement
    Parameter& operator=(Parameter&&);
    // \brief copy assignement
    Parameter& operator=(const Parameter&);
    // inheriting assignement operators
    using ParameterVariant::operator=;
    /*!
     * \brief constructor from C-string
     * \param[in] src: source
     */
    Parameter& operator=(const char* const);
    /*!
     * \brief constructor from a string_view
     * \param[in] src: source
     */
    Parameter& operator=(std::string_view);
    //
    ParameterVariant& as_std_variant() noexcept;
    const ParameterVariant& as_std_variant() const noexcept;
    //! \brief destructor
    ~Parameter();
  };  // end of struct Parameter

  //! \brief a simple alias
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
   * \return true if the given parameter has the given type
   * \param[in] p: parameters
   */
  template <typename ResultType>
  [[nodiscard]] bool is(const Parameter&) noexcept;

  //! \brief partial specialisation of the `is` function for double
  template <>
  [[nodiscard]] bool is<double>(const Parameter&) noexcept;
  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \throws if the parameter does not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] OptionalGetResultType<ResultType> get(
      Context&, const Parameter&) noexcept;

  //! \brief partial specialisation of the `get` function for double
  template <>
  [[nodiscard]] OptionalGetResultType<double> get<double>(
      Context&, const Parameter&) noexcept;
  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \throws if the parameter does not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing,
                                              const Parameter&);

  //! \brief partial specialisation of the `get` function for double
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing,
                                                  const Parameter&);

  /*!
   * \return true if the given parameter exists
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool contains(const Parameters&,
                                               std::string_view) noexcept;

  /*!
   * \return true if the given parameter has the given type
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists
   */
  template <typename ResultType>
  [[nodiscard]] bool is(attributes::Throwing,
                        const Parameters&,
                        std::string_view);
  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] OptionalGetResultType<ResultType> get(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \return value of the parameter
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] OptionalReference<const Parameter> get(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \return value of the parameter
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does not have the good type.
   */
  template <>
  [[nodiscard]] OptionalGetResultType<double> get<double>(
      Context&, const Parameters&, std::string_view) noexcept;
  /*!
   * \return value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \throws if the parameter exists but does not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] std::optional<ResultType> get_if(Context&,
                                                 const Parameters&,
                                                 std::string_view,
                                                 const ResultType&) noexcept;
  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in] throwing: dummy attribute to indicate
   * that this function may raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does
   * not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] GetResultType<ResultType> get(attributes::Throwing,
                                              const Parameters&,
                                              std::string_view);
  /*!
   * \return value of the parameter
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get(attributes::Throwing,
                                               const Parameters&,
                                               std::string_view);

  /*!
   * \return value of the parameter
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does not have the good type.
   */
  template <>
  [[nodiscard]] GetResultType<double> get<double>(attributes::Throwing,
                                                  const Parameters&,
                                                  std::string_view);

  /*!
   * \return value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * raise an exception
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \throws if the parameter exists but does not have the good type.
   */
  template <typename ResultType>
  [[nodiscard]] ResultType get_if(attributes::Throwing,
                                  const Parameters&,
                                  std::string_view,
                                  const ResultType&);
  /*!
   * \return value of the parameter if present, a default value otherwise
   *
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   */
  MFEM_MGIS_EXPORT [[nodiscard]] Parameter get_if(const Parameters&,
                                                  std::string_view,
                                                  const Parameter&) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.ixx"

#endif /* LIB_MFEM_MGIS_PARAMETER_HXX */
