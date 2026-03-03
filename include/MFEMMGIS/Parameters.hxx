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
   * \brief a structure representing a map associating a name to a parameter
   * \note inheritance is required here to allow forward declaration of the
   * Parameters structure.
   */
  struct MFEM_MGIS_EXPORT [[nodiscard]] Parameters
      : private std::map<std::string, Parameter, std::less<>> {
    /*!
     * \param[in, out] ctx: execution context
     * \param[in] k: key
     */
    static InvalidResult reportMissingKey(Context&, std::string_view) noexcept;
    /*!
     * \param[in, out] ctx: execution context
     * \param[in] n: parameter' name
     */
    static InvalidResult reportUnmatchedParameterType(
        Context&, std::string_view) noexcept;
    /*!
     * \brief throw an exception if the parameter type is not the expected one.
     * \param[in] n: name of the parameter
     */
    [[noreturn]] static void raiseUnmatchedParameterType(attributes::Throwing,
                                                         std::string_view);
    // exposing base class iterator
    using const_iterator =
        std::map<std::string, Parameter, std::less<>>::const_iterator;
    // inheriting constructors
    using std::map<std::string, Parameter, std::less<>>::map;
    //! \brief default constructor
    Parameters() noexcept;
    //! \brief copy constructor
    Parameters(const Parameters&) noexcept;
    //! \brief move constructor
    Parameters(Parameters&&) noexcept;
    //! \brief standard assignement
    Parameters& operator=(const Parameters&) noexcept;
    //! \brief move assignement
    Parameters& operator=(Parameters&&) noexcept;
    //! \return an iterator to the first element
    const_iterator begin() const noexcept;
    //! \return an iterator to the first element
    const_iterator cbegin() const noexcept;
    //! \return an iterator past the last element
    const_iterator end() const noexcept;
    //! \return an iterator past the last element
    const_iterator cend() const noexcept;
    /*!
     * \return true if the given parameter exists
     * \param[in] n: name
     */
    bool contains(std::string_view) const noexcept;
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(
        attributes::Throwing,
        const std::initializer_list<
            std::map<std::string, Parameter, std::less<>>::value_type>&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(attributes::Throwing, const Parameters&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(attributes::Throwing,
                       const std::map<std::string, Parameter>&);
    /*!
     * \brief insert a parameter using the given name
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     * \throws if the given parameter already exists
     */
    Parameters& insert(attributes::Throwing,
                       std::string_view,
                       const Parameter&);
    /*!
     * \brief insert a parameter using the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     */
    [[nodiscard]] bool insert(Context&,
                              std::string_view,
                              const Parameter&) noexcept;
    /*!
     * \return the parameter associated with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the parameter
     */
    OptionalReference<const Parameter> get(Context&,
                                           std::string_view) const noexcept;
    /*!
     * \return the parameter associated with the given name
     * \param[in] n: name of the parameter
     */
    const Parameter& get(attributes::Throwing, std::string_view) const;
    //! \brief destructor
    ~Parameters();
  };  // end of struct Parameters

  /*!
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] names: list of valid parameters names
   */
  MFEM_MGIS_EXPORT bool checkParameters(
      Context&, const Parameters&, const std::vector<std::string>&) noexcept;
  /*!
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] descriptions: descriptions of the allowed parameters
   */
  MFEM_MGIS_EXPORT bool checkParameters(
      Context&,
      const Parameters&,
      const std::map<std::string, std::string>&) noexcept;

  /*!
   * \param[in] parameters: parameters
   * \param[in] names: list of valid parameters names
   * \throws if an invalid parameter is present
   */
  MFEM_MGIS_EXPORT void checkParameters(attributes::Throwing,
                                        const Parameters&,
                                        const std::vector<std::string>&);
  /*!
   * \param[in] parameters: parameters
   * \param[in] descriptions: descriptions of the allowed parameters
   * \throws if an invalid parameter is present
   */
  MFEM_MGIS_EXPORT void checkParameters(
      attributes::Throwing,
      const Parameters&,
      const std::map<std::string, std::string>&);
  /*!
   * \brief extract the given parameters if they exists
   * \param[in, out] ctx: execution context
   * \param[in] parameters: parameters
   * \param[in] names: list of parameters names
   */
  MFEM_MGIS_EXPORT std::optional<Parameters> extract(
      Context&, const Parameters&, const std::vector<std::string>&) noexcept;
  /*!
   * \brief extract the given parameters if they exists
   * \param[in] parameters: parameters
   * \param[in] names: list of parameters names
   */
  MFEM_MGIS_EXPORT Parameters extract(attributes::Throwing,
                                      const Parameters&,
                                      const std::vector<std::string>&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameter.hxx"

#endif /* LIB_MFEM_MGIS_PARAMETERS_HXX */
