/*!
 * \file
 * \author Thomas Helfer
 * \date   07/12/2021
 * \brief header file declaring the
 * OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator class
 */

#ifndef LIB_MFEM_MGIS_ORTHOTROPICBIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ORTHOTROPICBIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX

#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief class implementing a behaviour integrator dedicated to
   * orthotropic micromorphic damage in two dimensions (the modelling
   * hypothesis has no effect on this specific behaviour as out of plane
   * damage gradients are assumed to be zero).
   */
  struct MFEM_MGIS_EXPORT
      OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator
      : BehaviourIntegratorBase {
    //! \brief a simple alias
    using RotationMatrix = std::array<real, 9u>;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator(
        const FiniteElementDiscretization &fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr);
    //
    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const override;
    real getIntegrationPointWeight(
        mfem::ElementTransformation &tr,
        const mfem::IntegrationPoint &ip) const noexcept override;
    /*!
     * \brief integrate the behaviour over the time step at the integration
     * points of the given element
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: nodal values of the micromorphic damage
     * \param[in] it: integration type
     * \return true on success
     */
    bool integrate(const mfem::FiniteElement &e,
                   mfem::ElementTransformation &tr,
                   const mfem::Vector &u,
                   const IntegrationType it) override;

    void updateResidual(mfem::Vector &Fe,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;

    /*!
     * \brief compute the contribution of the given element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns, unused
     * \note the stored tangent operator blocks are rotated in the global
     * frame in place.
     */
    void updateJacobian(mfem::DenseMatrix &Ke,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;

    void computeInnerForces(mfem::Vector &Fe,
                            const mfem::FiniteElement &e,
                            mfem::ElementTransformation &tr) override;
    //! \brief destructor
    ~OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator() override;

   private:
    /*!
     * \brief select the integration rule for the given element and element
     * transformation
     * \param[in] e: element
     * \param[in] t: element transformation
     * \return the integration rule
     */
    static const mfem::IntegrationRule &selectIntegrationRule(
        const mfem::FiniteElement &e, const mfem::ElementTransformation &t);
    /*!
     * \brief build the quadrature space for the given material
     * \param[in] throwing: throwing attribute.
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \return the partial quadrature space
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        attributes::Throwing throwing,
        const FiniteElementDiscretization &fed,
        const size_type m);

#ifndef MFEM_THREAD_SAFE
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif
  };  // end of struct
      // OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ORTHOTROPICBIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX \
        */
