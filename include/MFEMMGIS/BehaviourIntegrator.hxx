/*!
 * \file   MFEMMGIS/BehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  /*!
   * \brief abstract class for all behaviour integrators
   *
   * This class provides methods to:
   *
   * - integrate the behaviour over the time step
   * - compute the nodal forces du to the material reaction (see the
   *   `computeInnerForces` method).
   * - compute the stiffness matrix (see the `computeStiffnessMatrix` method).
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegrator {
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real) = 0;
    /*!
     * \brief compute the inner forces for the given element
     * \param[out] Fe: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     */
    virtual void computeInnerForces(mfem::Vector &,
                                    const mfem::FiniteElement &,
                                    mfem::ElementTransformation &,
                                    const mfem::Vector &) = 0;
    /*!
     * \brief compute the stiffness matrix for the given element
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     */
    virtual void computeStiffnessMatrix(mfem::DenseMatrix &,
                                        const mfem::FiniteElement &,
                                        mfem::ElementTransformation &,
                                        const mfem::Vector &) = 0;
    //! \brief destructor
    virtual ~BehaviourIntegrator();
  };  // end of struct BehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX */
