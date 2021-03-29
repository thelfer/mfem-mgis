/*!
 * \file   include/MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX
#define LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX

#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief a base class of the `StandardBehaviourIntegratorCRTPBase`
   * class to factorize code by static using the CRTP idiom.
   *
   * This class provides a way to optimise dynamic memory allocations.
   *
   * The `Child` class must provide:
   *
   * - a static method called `getIntegrationRule`
   * - a method called `updateGradients`
   * - a method called `updateInnerForces`
   * - a method called `updateStiffnessMatrix`
   * - a method called `getRotationMatrix`
   * - a method called `rotateGradients`
   * - a method called `rotateThermodynamicForces`
   * - a method called `rotateTangentOperatorBlocks`
   */
  template <typename Child>
  struct StandardBehaviourIntegratorCRTPBase : BehaviourIntegratorBase {
   protected:
    // inheriting `BehaviourIntegratorBase`' constructor
    using BehaviourIntegratorBase::BehaviourIntegratorBase;
    /*!
     * \brief compute the contribution of the element to the residual
     * \param[out] Fe: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     *
     * \note Thanks to the CRTP idiom, this implementation can call
     * the `updateGradients` and the `updateInnerForces` methods defined
     * in the derived class without a virtual call. Those call may
     * even be inlined.
     * \note The implementation of the `updateResidual` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the `implementUpdateResidual` is
     * generated.
     */
    void implementUpdateResidual(mfem::Vector &,
                                 const mfem::FiniteElement &,
                                 mfem::ElementTransformation &,
                                 const mfem::Vector &);
    /*!
     * \brief compute the contribution of the element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     *
     * \note The implementation of the `updateJacobian` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the
     * `implementUpdateJacobian` is generated.
     */
    void implementUpdateJacobian(mfem::DenseMatrix &,
                                 const mfem::FiniteElement &,
                                 mfem::ElementTransformation &);
    /*!
     * \brief compute the contribution of the element to the inner forces
     * \param[out] Fe: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     */
    void implementComputeInnerForces(mfem::Vector&,
                                     const mfem::FiniteElement &,
                                     mfem::ElementTransformation &);
    //! \brief destructor
    ~StandardBehaviourIntegratorCRTPBase() override;

#ifndef MFEM_THREAD_SAFE
   private:
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif

  };  // end of StandardBehaviourIntegratorCRTPBase

}  // end of namespace mfem_mgis

#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.ixx"

#endif /* LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX */
