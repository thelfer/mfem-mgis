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
                                        Parameters,
                                        std::function<real(const real)>>;

  /*!
   *
   */
  struct MFEM_MGIS_EXPORT Parameter : ParameterVariant {
    /*!
     * \brief throw an exception if the parameter type is not the expected one.
     */
    [[noreturn]] static void raiseUnmatchedParameterType();
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
    //! \brief destructor
    ~Parameter();
  };  // end of struct Parameter

  //! \brief a simple alias
  template <typename ResultType>
  using GetResultType = std::conditional_t<std::is_same_v<ResultType, double>,
                                           ResultType,
                                           const ResultType&>;

  /*!
   * \return true if the given parameter has the given type
   * \param[in] p: parameters
   */
  template <typename ResultType>
  bool is(const Parameter&);

  //! \brief partial specialisation of the `is` function for double
  template <>
  bool is<double>(const Parameter&);

  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] p: parameters
   * \throws if the parameter does not have the good type.
   */
  template <typename ResultType>
  GetResultType<ResultType> get(const Parameter&);

  //! \brief partial specialisation of the `get` function for double
  template <>
  GetResultType<double> get<double>(const Parameter&);

  /*!
   * \return true if the given parameter exists
   * \param[in] p: parameters
   * \param[in] n: name
   */
  bool contains(const Parameters&, std::string_view);

  /*!
   * \return true if the given parameter has the given type
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists
   */
  template <typename ResultType>
  bool is(const Parameters&, std::string_view);

  /*!
   * \return value of the parameter
   * \tparam ResultType: expected type of the parameter
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does not have the good type.
   */
  template <typename ResultType>
  GetResultType<ResultType> get(const Parameters&, std::string_view);

  /*!
   * \return value of the parameter
   * \param[in] p: parameters
   * \param[in] n: name
   * \throws if the parameter does not exists or does not have the good type.
   */
  template <>
  GetResultType<double> get<double>(const Parameters&, std::string_view);

  /*!
   * \return value of the parameter if present, a default value otherwise
   * \tparam ResultType: expected type of the parameter
   * \param[in] p: parameters
   * \param[in] n: name
   * \param[in] v: default value
   * \throws if the parameter exists but does not have the good type.
   */
  template <typename ResultType>
  ResultType get_if(const Parameters&, std::string_view, const ResultType&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.ixx"

#endif /* LIB_MFEM_MGIS_PARAMETER_HXX */
