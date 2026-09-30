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
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  [[nodiscard]] MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeVonMisesEquivalentStress(Context& ctx,
                                  const Material& m,
                                  const Material::StateSelection s);
  /*!
   * \brief compute the von Mises equivalent stress.
   *
   * \note For finite strain behaviours, the von Mises equivalent stress of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] seq: von Mises equivalent stress
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  [[nodiscard]] MFEM_MGIS_EXPORT bool computeVonMisesEquivalentStress(
      Context& ctx,
      PartialQuadratureFunction& seq,
      const Material& m,
      const Material::StateSelection s);

  /*!
   * \brief compute the eigen values of the stress.
   *
   * \note For finite strain behaviours, the eigen values of the
   * Cauchy stress are returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return the eigen values of the stress
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeEigenStresses(Context& ctx,
                       const Material& m,
                       const Material::StateSelection s);
  /*!
   * \brief compute the eigen values of the stress.
   *
   * \note For finite strain behaviours, the eigen values of the
   * Cauchy stress are returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] svp: eigen values of the stress
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT bool computeEigenStresses(Context& ctx,
                                             PartialQuadratureFunction& svp,
                                             const Material& m,
                                             const Material::StateSelection s);

  /*!
   * \brief compute the first (maximum) eigen value of the stress.
   *
   * \note For finite strain behaviours, the first eigen value of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return the first eigen value of the stress
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeFirstEigenStress(Context& ctx,
                          const Material& m,
                          const Material::StateSelection s);
  /*!
   * \brief compute the first (maximum) eigen value of the stress.
   *
   * \note For finite strain behaviours, the first eigen value of the
   * Cauchy stress is returned.
   * \note This function currently does not work in plane stress for finite
   * strain behaviours.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] s1: first eigen value of the stress
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT bool computeFirstEigenStress(
      Context& ctx,
      PartialQuadratureFunction& s1,
      const Material& m,
      const Material::StateSelection s);
  /*!
   * \brief compute the stress in the global frame.
   *
   * \note For finite strain behaviours, the first Piola-Kirchhoff stress
   * tensor is returned.
   * \note This function is only valid for orthotropic behaviours.
   *
   * \return the stress in the global frame
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeStressInGlobalFrame(Context& ctx,
                             const Material& m,
                             const Material::StateSelection s);
  /*!
   * \brief compute the stress in the global frame.
   *
   * \note For finite strain behaviours, the first Piola-Kirchhoff stress
   * tensor is returned.
   * \note This function is only valid for orthotropic behaviours.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] rstress: stress in the global frame
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT bool computeStressInGlobalFrame(
      Context& ctx,
      PartialQuadratureFunction& rstress,
      const Material& m,
      const Material::StateSelection s);

  /*!
   * \brief compute the Cauchy stress in the global frame.
   *
   * \note This function is only valid for finite strain behaviours
   * \note This function currently does not work in plane stress.
   * \note This function is only valid for orthotropic behaviours.
   *
   * \return the Cauchy stress in the global frame
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction>
  computeCauchyStressInGlobalFrame(Context& ctx,
                                   const Material& m,
                                   const Material::StateSelection s);
  /*!
   * \brief compute the Cauchy stress in the global frame.
   *
   * \note This function is only valid for finite strain behaviours
   * \note This function currently does not work in plane stress.
   * \note This function is only valid for orthotropic behaviours.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] sig: Cauchy stress in the global frame
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT bool computeCauchyStressInGlobalFrame(
      Context& ctx,
      PartialQuadratureFunction& sig,
      const Material& m,
      const Material::StateSelection s);

  /*!
   * \brief compute the Cauchy stress.
   *
   * \note This function is only valid for finite strain behaviours
   * \note This function currently does not work in plane stress.
   *
   * \return the Cauchy stress
   * \param[in, out] ctx: execution context
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction> computeCauchyStress(
      Context& ctx, const Material& m, const Material::StateSelection s);
  /*!
   * \brief compute the Cauchy stress.
   *
   * \note This function is only valid for finite strain behaviours
   * \note This function currently does not work in plane stress.
   *
   * \return true on success
   * \param[in, out] ctx: execution context
   * \param[out] sig: Cauchy stress
   * \param[in] m: material
   * \param[in] s: selection of the state considered (beginning of time step,
   * end of time step)
   */
  MFEM_MGIS_EXPORT bool computeCauchyStress(Context& ctx,
                                            PartialQuadratureFunction& sig,
                                            const Material& m,
                                            const Material::StateSelection s);

#endif /* MGIS_FUNCTION_SUPPORT */

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_MECHANICALPOSTPROCESSINGS_HXX */
