/*!
 * \file   include/MFEMMGIS/Material.ixx
 * \brief
 * \author Thomas Helfer
 * \date   03/03/2021
 */

#ifndef LIB_MFEM_MGIS_MATERIAL_IXX
#define LIB_MFEM_MGIS_MATERIAL_IXX

namespace mfem_mgis {

  inline std::array<real, 9u> Material::getRotationMatrixAtIntegrationPoint(
      const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrixAtIntegrationPoint

  inline mgis::behaviour::MaterialStateManager &getStateManager(
      Material &m, const Material::StateSelection s) {
    if (s == Material::END_OF_TIME_STEP) {
      return m.s1;
    }
    return m.s0;
  }

  inline const mgis::behaviour::MaterialStateManager &getStateManager(
      const Material &m, const Material::StateSelection s) {
    if (s == Material::END_OF_TIME_STEP) {
      return m.s1;
    }
    return m.s0;
  }

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATERIAL_IXX */
