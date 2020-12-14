/*!
 * \file   src/SmallStrainMechanicalBehaviourIntegrator.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.hxx"

namespace mfem_mgis {

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

}  // end of namespace mfem_mgis



