/*!
 * \file MFEMMGIS/FBarBehaviourIntegrators.hxx
 * \brief
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

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFBarBehaviourIntegrators(
      Context &,
      const FiniteElementDiscretization &,
      const size_type,
      std::unique_ptr<const Behaviour>,
      const Parameters &) noexcept;

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressFBarBehaviourIntegrators(
      Context &,
      const FiniteElementDiscretization &,
      const size_type,
      std::unique_ptr<const Behaviour>,
      const Parameters &) noexcept;

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFBarBehaviourIntegrators(
      Context &,
      const FiniteElementDiscretization &,
      const size_type,
      std::unique_ptr<const Behaviour>,
      const Parameters &) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORS_HXX */
