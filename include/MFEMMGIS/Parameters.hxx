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
#include "MFEMMGIS/Parameter.hxx"

namespace mfem_mgis {

  /*!
   * \brief a structure representing a map associating a name to a parameter
   * \note inheritance is required here to allow forward declaration of the
   * Parameters structure.
   */
  struct Parameters : private std::map<std::string, Parameter, std::less<>> {
    /*!
     * \brief throw an exception if the parameter type is not the expected one.
     * \param[in] p: parameter
     */
    [[noreturn]] static void raiseUnmatchedParameterType(std::string_view);
    // exposing base class iterator
    using const_iterator =
        std::map<std::string, Parameter, std::less<>>::const_iterator;
    // inheriting constructors
    using std::map<std::string, Parameter, std::less<>>::map;
    //! \brief default constructor
    Parameters();
    //! \brief copy constructor
    Parameters(const Parameters&);
    //! \brief move constructor
    Parameters(Parameters&&);
    //! \brief standard assignement
    Parameters& operator=(const Parameters&);
    //! \brief move assignement
    Parameters& operator=(Parameters&&);
    //! \return an iterator to the first element
    const_iterator begin() const;
    //! \return an iterator to the first element
    const_iterator cbegin() const;
    //! \return an iterator past the last element
    const_iterator end() const;
    //! \return an iterator past the last element
    const_iterator cend() const;
    /*!
     * \return true if the given parameter exists
     * \param[in] n: name
     */
    bool contains(std::string_view) const;
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(
        const std::initializer_list<
            std::map<std::string, Parameter, std::less<>>::value_type>&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(const Parameters&);
    /*!
     * \brief insert parameters
     * \param[in] src: parameters
     * \throws if one of the parameters already exists
     */
    Parameters& insert(const std::map<std::string, Parameter>&);
    /*!
     * \brief insert a parameter using the given name
     * \param[in] n: name of the parameter
     * \param[in] p: parameter
     * \throws if the given parameter already exists
     */
    Parameters& insert(std::string_view, const Parameter&);
    /*!
     * \return the parameter associated with the given name
     * \param[in] n: name of the parameter
     */
    const Parameter& get(std::string_view) const;
    //! \brief destructor
    ~Parameters();
  };  // end of struct Parameters

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
  const ResultType& get(const Parameters&, std::string_view);

  /*!
   * \param[in] p: parameters
   * \param[in] names: list of valid parameters names
   * \throws if an invalid parameter is present
   */
  void checkParameters(const Parameters&, const std::vector<std::string>&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Parameters.ixx"

#endif LIB_MFEM_MGIS_PARAMETERS_HXX
