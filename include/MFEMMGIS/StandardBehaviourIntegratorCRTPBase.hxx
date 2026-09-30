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
   * \brief a base class for standard behaviour integrators, factorizing code
   * by static polymorphism using the CRTP idiom.
   *
   * This class provides a way to optimise dynamic memory allocations.
   *
   * The `Child` class must provide:
   *
   * - a method called `getIntegrationRule`
   * - a method called `getIntegrationPointWeight`
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
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    StandardBehaviourIntegratorCRTPBase(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::unique_ptr<const Behaviour> b_ptr);
    /*!
     * \brief integrate the behaviour over the time step
     *
     * If successful, the values of the thermodynamic forces, consistent
     * tangent operator and internal state variables are updated.
     *
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \return true on success
     */
    bool implementIntegrate(const mfem::FiniteElement& e,
                            mfem::ElementTransformation& tr,
                            const mfem::Vector& u,
                            const IntegrationType it);
    /*!
     * \brief compute the contribution of the element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns, unused
     *
     * \note Thanks to the CRTP idiom, this implementation can call
     * the `updateInnerForces` method defined
     * in the derived class without a virtual call. This call may
     * even be inlined.
     * \note The implementation of the `updateResidual` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the `implementUpdateResidual` is
     * generated.
     */
    void implementUpdateResidual(mfem::Vector& Fe,
                                 const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr,
                                 const mfem::Vector& u);
    /*!
     * \brief compute the contribution of the element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     *
     * \note The implementation of the `updateJacobian` in the
     * `Child` class trivially calls this method. This indirection is made to
     * control where the code associated to the
     * `implementUpdateJacobian` is generated.
     */
    void implementUpdateJacobian(mfem::DenseMatrix& Ke,
                                 const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr);
    /*!
     * \brief compute the contribution of the element to the inner forces
     * \param[out] Fe: inner forces
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     */
    void implementComputeInnerForces(mfem::Vector& Fe,
                                     const mfem::FiniteElement& e,
                                     mfem::ElementTransformation& tr);
    //! \brief destructor
    ~StandardBehaviourIntegratorCRTPBase() override;

#ifndef MFEM_THREAD_SAFE
   protected:
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif

  };  // end of StandardBehaviourIntegratorCRTPBase

}  // end of namespace mfem_mgis

#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.ixx"

#endif /* LIB_MFEM_MGIS_ISOTROPICSTANDARDBEHAVIOURINTEGRATORCRTPBASE_HXX */
