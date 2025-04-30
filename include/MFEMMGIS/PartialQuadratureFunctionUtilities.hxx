/*!
 * \file   MFEMMIGS/PartialQuadratureFunctionUtilities.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   30/04/2025
 */

#ifndef LIB_MFEMMIGS_PARTIALQUADRATUREFUNCTIONUTILITIES_HXX
#define LIB_MFEMMIGS_PARTIALQUADRATUREFUNCTIONUTILITIES_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis{

  /*!
   * \brief rotate the thermodynamic forces in the global frame
   * \param[out] f: quadrature function containing the thermodynamic forces in the global frame
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool rotateThermodynamicsForcesToGlobalFrame(
      PartialQuadratureFunction &,
      Material &,
      const Material::StateSelection = Material::END_OF_TIME_STEP);

} // end of namespace mfem_mgis

#endif /* LIB_MFEMMIGS_PARTIALQUADRATUREFUNCTIONUTILITIES_HXX */
