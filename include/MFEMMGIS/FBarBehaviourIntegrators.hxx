/*!
 * \file MFEMMGIS/FBarBehaviourIntegrators.hxx
 * \brief  This file declares the functions generating the FBar behaviour
 * integrators
 * \author Thomas Helfer
 * \date   22/03/2026
 */

#ifndef LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORS_HXX
#define LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORS_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief generate an FBar behaviour integrator under the plane strain
   * hypothesis
   * \param[in, out] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] m: material attribute
   * \param[in] b: behaviour
   * \param[in] params: additional parameters, none expected
   * \return the new behaviour integrator, a null pointer on failure
   */
  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFBarBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

  /*!
   * \brief generate an FBar behaviour integrator under the tridimensional
   * hypothesis
   * \param[in, out] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] m: material attribute
   * \param[in] b: behaviour
   * \param[in] params: additional parameters, none expected
   * \return the new behaviour integrator, a null pointer on failure
   */
  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFBarBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORS_HXX */
