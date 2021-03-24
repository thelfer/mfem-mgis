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

namespace mfem_mgis {

  //! \brief a simple alias
  using Parameter =
      std::variant<std::monostate, int, real, std::function<real(const real)>>;

}  // end of namespace mfem_mgis

#endif LIB_MFEM_MGIS_PARAMETER_HXX
