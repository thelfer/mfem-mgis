/*!
 * \file   MFEMMGIS/MechanicalPostProcessings.hxx
 * \brief
 * \author Thomas Helfer
 * \date   22/05/2025
 */

#ifndef LIB_MFEMMGIS_MECHANICALPOSTPROCESSINGS_HXX
#define LIB_MFEMMGIS_MECHANICALPOSTPROCESSINGS_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

#ifdef MGIS_FUNCTION_SUPPORT

  /*!
   * \brief compute the von Mises equivalent stress.
   *
   * \note For finite strain behaviours, the von Mises equivalent stress of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return the von Mises stress
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginnig of time step, end
   * of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStress(Context&,
                                  const Material&,
                                  const Material::StateSelection);

  /*!
   * \brief compute the eigen values of the stress.
   *
   * \note For finite strain behaviours, the eigen values of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return the eigen values of the stress
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginnig of time step, end
   * of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeEigenStresses(Context&,
                       const Material&,
                       const Material::StateSelection);

  /*!
   * \brief compute the first (maximum) eigen value of the stress.
   *
   * \note For finite strain behaviours, the first eigen value of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return the eigen values of the stress
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginnig of time step, end
   * of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeFirstEigenStress(Context&,
                          const Material&,
                          const Material::StateSelection);

  /*!
   * \brief compute the Cauchy stress in the global frame.
   *
   * \note This function is only valid for finite strain behaviours
   * \note This function currently does not work in plane stress.
   *
   * \return the eigen values of the stress
   * \param[in] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginnig of time step, end
   * of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeCauchyStressInGlobalFrame(Context&,
                                   const Material&,
                                   const Material::StateSelection);

#endif /* MGIS_FUNCTION_SUPPORT */

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_MECHANICALPOSTPROCESSINGS_HXX */
