/*!
 * \file   SmallStrainMechanicalBehaviourIntegrator.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#ifndef LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX
#define LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX

#include <utility>
#include "MGIS/Raise.hxx"

namespace mfem_mgis {

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<H>::
      SmallStrainMechanicalBehaviourIntegrator(const mfem::FiniteElementSpace & fs,
                                               const size_type m,
                                               std::unique_ptr<const Behaviour> b_ptr)
      : SmallStrainMechanicalBehaviourIntegratorBase(fs, m, std::move(b_ptr)) {
    this->checkHypotheses(H);
  }  // end of SmallStrainMechanicalBehaviourIntegrator

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::computeInnerForces(
      const mfem::FiniteElement &,
      mfem::ElementTransformation &,
      const mfem::Vector &,
      mfem::Vector &) {
  }  // end of computeInnerForces

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::computeStiffnessMatrix(
      const mfem::FiniteElement &,
      mfem::ElementTransformation &,
      const mfem::Vector &,
      mfem::DenseMatrix &) {}  // end of computeStiffnessMatrix

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<
      H>::~SmallStrainMechanicalBehaviourIntegrator() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX */
