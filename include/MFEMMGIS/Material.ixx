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
      Material &m, const Material::StateSelection s) noexcept {
    if (s == Material::END_OF_TIME_STEP) {
      return m.s1;
    }
    return m.s0;
  }

  inline const mgis::behaviour::MaterialStateManager &getStateManager(
      const Material &m, const Material::StateSelection s) noexcept {
    if (s == Material::END_OF_TIME_STEP) {
      return m.s1;
    }
    return m.s0;
  }

#ifdef MGIS_FUNCTION_SUPPORT

  inline bool check(AbstractErrorHandler &eh,
                    const RotationMatrixEvaluator &e) {
    return e.check(eh);
  }  // end of check

  constexpr mgis::size_type getNumberOfComponents(
      const RotationMatrixEvaluator &) noexcept {
    return 9u;
  }  // end of getNumberOfComponents

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATERIAL_IXX */
