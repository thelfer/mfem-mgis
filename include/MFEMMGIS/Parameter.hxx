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
#include <string_view>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  //! \brief a simple alias
  using ParameterVariant = std::variant<std::monostate,
                                        bool,
                                        int,
                                        real,
                                        std::string,
                                        std::function<real(const real)>>;

  /*!
   *
   */
  struct MFEM_MGIS_EXPORT Parameter : ParameterVariant {
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

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PARAMETER_HXX */
