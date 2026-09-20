/*!
 * \file   include/MFEMMGIS/Parameters.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEM_MGIS_PARAMETERS_HXX
#define LIB_MFEM_MGIS_PARAMETERS_HXX

#include <map>
#include <vector>
#include <string>
#include <string_view>
#include <initializer_list>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;

  /*!
   * \brief A structure representing a map associating a name to a parameter.
   * \note Inheritance is required here to allow forward declaration of the
   * Parameters structure.
   */
  struct MFEM_MGIS_EXPORT [[nodiscard]] Parameters
      : private std::map<std::string, Parameter, std::less<>> {
    /*!
     * \brief Report that a key is missing.
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     */
    static InvalidResult reportMissingKey(Context&, std::string_view) noexcept;
    /*!
     * \brief Report that the type of a parameter is not the expected one.
     * \param[in, out] ctx: execution context
     * \param[in] n: parameter' name
     */
    static InvalidResult reportUnmatchedParameterType(
        Context&, std::string_view) noexcept;
    /*!
     * \brief Throw an exception if the parameter type is not the expected one.
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     * \param[in] n: name of the parameter
     */
    [[noreturn]] static void raiseUnmatchedParameterType(attributes::Throwing,
                                                         std::string_view);
    // exposing base class iterator
    using const_iterator =
        std::map<std::string, Parameter, std::less<>>::const_iterator;
    // inheriting constructors
    using std::map<std::string, Parameter, std::less<>>::map;
    //! \brief Default constructor
    Parameters() noexcept;
    //! \brief Copy constructor
    Parameters(const Parameters&) noexcept;
    //! \brief Move constructor
    Parameters(Parameters&&) noexcept;
    //! \brief Standard assignment
    Parameters& operator=(const Parameters&) noexcept;
    //! \brief Move assignment
    Parameters& operator=(Parameters&&) noexcept;
    //! \return An iterator to the first element
    const_iterator begin() const noexcept;
    //! \return An iterator to the first element
    const_iterator cbegin() const noexcept;
    //! \return An iterator past the last element
    const_iterator end() const noexcept;
    //! \return An iterator past the last element
    const_iterator cend() const noexcept;
    /*!
     * \return True if the given parameter exists
     * \param[in] n: name
     */
    bool contains(std::string_view) const noexcept;
    /*!
     * \brief Insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     */
    Parameters& insert(
        attributes::Throwing,
        const std::initializer_list<
            std::map<std::string, Parameter, std::less<>>::value_type>&);
    /*!
     * \brief Insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     */
    Parameters& insert(attributes::Throwing, const Parameters&);
    /*!
     * \brief Insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     */
    Parameters& insert(attributes::Throwing,
                       const std::map<std::string, Parameter>&);
    /*!
     * \brief Insert a parameter using the given name
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     * \throws if the given parameter already exists
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     */
    Parameters& insert(attributes::Throwing,
                       std::string_view,
                       const Parameter&);
    /*!
     * \brief Insert a parameter using the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     */
    [[nodiscard]] bool insert(Context&,
                              std::string_view,
                              const Parameter&) noexcept;
    /*!
     * \return The parameter associated with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the parameter
     */
    OptionalReference<const Parameter> get(Context&,
                                           std::string_view) const noexcept;
    /*!
     * \return The parameter associated with the given name
     * \param[in] n: name of the parameter
     * \throws if the parameter does not exist
     * \note this function shall only be used in constructors or functions with
     * the `attributes::Throwing` attribute when no context is available.
     */
    const Parameter& get(attributes::Throwing, std::string_view) const;
    /*!
     * \brief Replace the given parameter by the new value
     * \param[in] n: name of the parameter
     * \param[in] v: new value
     */
    Parameters& replaceOrInsert(std::string_view, const Parameter&) noexcept;
    //
    using std::map<std::string, Parameter, std::less<>>::size;
    using std::map<std::string, Parameter, std::less<>>::empty;
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    [[deprecated]] Parameters& insert(
        const std::initializer_list<
            std::map<std::string, Parameter, std::less<>>::value_type>&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    [[deprecated]] Parameters& insert(const Parameters&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    [[deprecated]] Parameters& insert(const std::map<std::string, Parameter>&);
    /*!
     * \brief insert a parameter using the given name
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     * \throws if the given parameter already exists
     */
    [[deprecated]] Parameters& insert(std::string_view, const Parameter&);
    //! \brief destructor
    ~Parameters();
  };  // end of struct Parameters

  /*!
   * \brief Check if the given parameters are valid
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] names: list of valid parameters names
   * \return True if all parameters are valid
   *
   * \note prefer using `ParametersValidator`
   */
  MFEM_MGIS_EXPORT bool checkParameters(
      Context&, const Parameters&, const std::vector<std::string>&) noexcept;
  /*!
   * \brief Check if the given parameters are valid
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] descriptions: descriptions of the allowed parameters
   * \return True if all parameters are valid
   *
   * \note prefer using `ParametersValidator`
   */
  MFEM_MGIS_EXPORT bool checkParameters(
      Context&,
      const Parameters&,
      const std::map<std::string, std::string>&) noexcept;

  /*!
   * \brief Check if the given parameters are valid
   * \param[in] parameters: parameters
   * \param[in] names: list of valid parameters names
   * \throws if an invalid parameter is present
   *
   * \note prefer using `ParametersValidator`
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT void checkParameters(attributes::Throwing,
                                        const Parameters&,
                                        const std::vector<std::string>&);
  /*!
   * \brief Check if the given parameters are valid
   * \param[in] parameters: parameters
   * \param[in] descriptions: descriptions of the allowed parameters
   * \throws if an invalid parameter is present
   *
   * \note prefer using `ParametersValidator`
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT void checkParameters(
      attributes::Throwing,
      const Parameters&,
      const std::map<std::string, std::string>&);
  /*!
   * \brief Extract the given parameters if they exist
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] names: list of parameters names
   * \return The extracted parameters if they exist
   */
  MFEM_MGIS_EXPORT std::optional<Parameters> extract(
      Context&, const Parameters&, const std::vector<std::string>&) noexcept;
  /*!
   * \brief Extract the given parameters if they exist
   * \param[in] parameters: parameters
   * \param[in] names: list of parameters names
   * \throws if an error occurs
   *
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT Parameters extract(attributes::Throwing,
                                      const Parameters&,
                                      const std::vector<std::string>&);
  /*!
   * \brief Extract the given parameters if they exist
   * \param[in] parameters: parameters
   * \param[in] descriptions: description of parameters to be extracted
   * \throws if an error occurs
   *
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT Parameters
  extract(attributes::Throwing,
          const Parameters&,
          const std::map<std::string, std::string>&);
  /*!
   * \return The information required to build an object from a factory
   * \param[in, out] ctx: execution context
   * \param[in] p: parameters
   */
  MFEM_MGIS_EXPORT std::optional<std::pair<std::string, Parameters>>
  extractFactoryArgument(Context& ctx, const Parameters& parameters) noexcept;
  /*!
   * \return The information required to build an object from a factory
   * \param[in] p: parameters
   * \throws if an error occurs
   *
   * \note this function shall only be used in constructors or functions with
   * the `attributes::Throwing` attribute when no context is available.
   */
  MFEM_MGIS_EXPORT std::pair<std::string, Parameters> extractFactoryArgument(
      attributes::Throwing, const Parameters&);
  /*!
   * \brief Remove in the output parameters the given parameters if they exist
   * \param[in] parameters: parameters
   * \param[in] names: list of parameters names to be removed
   * \return The parameters with the given names removed
   */
  MFEM_MGIS_EXPORT Parameters remove(const Parameters&,
                                     const std::vector<std::string>&) noexcept;
  /*!
   * \brief Remove in the output parameters the given parameters if they exist
   * \param[in] parameters: parameters
   * \param[in] descriptions: description of parameters to be extracted
   * \return The parameters with the given names removed
   */
  MFEM_MGIS_EXPORT Parameters
  remove(const Parameters&, const std::map<std::string, std::string>&) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.hxx"

#endif /* LIB_MFEM_MGIS_PARAMETERS_HXX */
