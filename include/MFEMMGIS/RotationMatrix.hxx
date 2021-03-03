/*!
 * \file   include/MFEMMGIS/RotationMatrix.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_ROTATIONMATRIX_HXX
#define LIB_MFEM_MGIS_ROTATIONMATRIX_HXX

#include <array>
#include <tuple>
#include <variant>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  /*!
   * \brief definition of the rotation matrix in 2D
   * The rotation matrix is defined by:
   * - a constant matrix, computed once for all
   * - the definition of the first direction of orthotropy
   */
  using RotationMatrix2D =
      std::variant<std::monostate,
                   std::array<real, 9u>,
                   std::shared_ptr<PartialQuadratureFunction>>;
  //! a simple alias
  using MaterialAxis3D =
      std::variant<std::array<real, 3u>,
                   std::shared_ptr<PartialQuadratureFunction>>;
  /*!
   * \brief definition of the rotation matrix in 3D
   * The rotation matrix is defined by:
   * - a constant matrix, computed once for all
   * - the two directions of orthotropy which can be given either as a
   * partial function or as a fixed direction.
   */
  using RotationMatrix3D = std::variant<std::monostate,
                                        std::array<real, 9u>,
                                        std::array<MaterialAxis3D, 2u>>;

};  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ROTATIONMATRIX_HXX */
