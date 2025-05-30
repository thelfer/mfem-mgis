/*!
 * \file include/MFEMMGIS/BidimensionalMicromorphicDamageBehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/12/2021
 * \brief header file declaring the
 * BidimensionalMicromorphicDamageBehaviourIntegrator class
 */

#ifndef LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX

#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief class implementing the a behaviour integrator dedicated to
   * micromorphic damage in two dimensions (the modelling hypothesis has no
   * effect on this specific behaviour as out of plane damage gradients are
   * assumed to be zero).
   */
  struct MFEM_MGIS_EXPORT BidimensionalMicromorphicDamageBehaviourIntegrator
      : BehaviourIntegratorBase {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    BidimensionalMicromorphicDamageBehaviourIntegrator(
        const FiniteElementDiscretization &,
        const size_type,
        std::unique_ptr<const Behaviour>);
    //
    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &,
        const mfem::ElementTransformation &) const override;
    real getIntegrationPointWeight(
        mfem::ElementTransformation &,
        const mfem::IntegrationPoint &) const noexcept override;
    bool integrate(const mfem::FiniteElement &,
                   mfem::ElementTransformation &,
                   const mfem::Vector &,
                   const IntegrationType) override;

    void updateResidual(mfem::Vector &,
                        const mfem::FiniteElement &,
                        mfem::ElementTransformation &,
                        const mfem::Vector &) override;

    void updateJacobian(mfem::DenseMatrix &,
                        const mfem::FiniteElement &,
                        mfem::ElementTransformation &,
                        const mfem::Vector &) override;

    void computeInnerForces(mfem::Vector &,
                            const mfem::FiniteElement &,
                            mfem::ElementTransformation &) override;
    //! \brief destructor
    ~BidimensionalMicromorphicDamageBehaviourIntegrator() override;

   private:
    /*!
     * \return the integration rule for the given element and  * element
     * transformation. \param[in] e: element \param[in] tr: element
     * transformation
     */
    static const mfem::IntegrationRule &selectIntegrationRule(
        const mfem::FiniteElement &, const mfem::ElementTransformation &);
    /*!
     * \brief build the quadrature space for the given  * material
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        const FiniteElementDiscretization &, const size_type);

#ifndef MFEM_THREAD_SAFE
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif
  };  // end of struct BidimensionalMicromorphicDamageBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BIDIMENSIONALMICROMORPHICDAMAGEBEHAVIOURINTEGRATOR_HXX \
        */
