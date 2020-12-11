/*!
 * \file   include/MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#ifndef LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX

#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for behaviour integrators based on small strain
   * mechanical behaviours.
   */
  struct SmallStrainMechanicalBehaviourIntegratorBase
      : BehaviourIntegratorBase {
    /*!
     * \brief constructor
     * \param[in] fs: finite element space.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    SmallStrainMechanicalBehaviourIntegratorBase(
        const mfem::FiniteElementSpace &,
        const size_type,
        std::shared_ptr<const Behaviour>);
    /*!
     * \return the integration rule for the given element and element
     * transformation.
     * \param[in] e: element
     * \param[in] tr: element transformation
     */
    virtual const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &, const mfem::ElementTransformation &) const;

    //! \brief destructor
    ~SmallStrainMechanicalBehaviourIntegratorBase() override;

  };  // end of struct SmallStrainMechanicalBehaviourIntegratorBase

  /*!
   * \brief a base class of the `SmallStrainMechanicalBehaviourIntegrator` class
   * to factorize code by static using the CRTP idiom.
   *
   * This class provides a way to optimise memory allocations
   */
  template <typename Child>
  struct SmallStrainMechanicalBehaviourIntegratorCRTPBase
      : public SmallStrainMechanicalBehaviourIntegratorBase {
    // inheriting constructors
    using SmallStrainMechanicalBehaviourIntegratorBase::
        SmallStrainMechanicalBehaviourIntegratorBase;
    //! \brief destructor
    ~SmallStrainMechanicalBehaviourIntegratorCRTPBase() override;

   protected:
    /*!
     * \brief compute the inner forces for the given element
     * \param[out] Fe: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     *
     * \note Thanks to the CRTP idiom, this implementation can call
     * the `updateStrain` and the `updateInnerForces` methods defined
     * in the derived class without a virtual call. Those call may
     * even be inlined.
     * \note The implementation of the `computeInnerForces` in the
     * `SmallStrainMechanicalBehaviourIntegrator` class trivially
     * calls this method. This indirection is made to control where
     * the code associated to the `implementComputeInnerForces`
     * is generated.
     */
    void implementComputeInnerForces(mfem::Vector &,
                                     const mfem::FiniteElement &,
                                     mfem::ElementTransformation &,
                                     const mfem::Vector &);
    /*!
     * \brief compute the stiffness matrix for the given element
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     *
     * \note This method is meant to implement the `computeStiffnessMatrix`
     * method. It is doned in
     */
    void implementComputeStiffnessMatrix(mfem::DenseMatrix &,
                                         const mfem::FiniteElement &,
                                         mfem::ElementTransformation &);

#ifndef MFEM_THREAD_SAFE
   private:
    //! \brief matrix used to store the derivatives of the shape functions
    mfem::DenseMatrix dshape;
#endif

  };  // end of SmallStrainMechanicalBehaviourIntegratorCRTPBase

  /*!
   * \brief behaviour integrators based on small strain
   * mechanical behaviours.
   *
   * \param[in] H: modelling hypothesis
   */
  template <Hypothesis H>
  struct MFEM_MGIS_EXPORT SmallStrainMechanicalBehaviourIntegrator
      : SmallStrainMechanicalBehaviourIntegratorCRTPBase<
            SmallStrainMechanicalBehaviourIntegrator<H>> {
    /*!
     * \brief constructor
     * \param[in] fs: finite element space.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    SmallStrainMechanicalBehaviourIntegrator(const mfem::FiniteElementSpace &,
                                             const size_type,
                                             std::shared_ptr<const Behaviour>);

    void computeInnerForces(mfem::Vector &,
                            const mfem::FiniteElement &,
                            mfem::ElementTransformation &,
                            const mfem::Vector &) override;

    void computeStiffnessMatrix(mfem::DenseMatrix &,
                                const mfem::FiniteElement &,
                                mfem::ElementTransformation &,
                                const mfem::Vector &) override;

    //! \brief destructor
    ~SmallStrainMechanicalBehaviourIntegrator() override;

   protected:
    //! \brief allow the CRTP base class the protected members
    friend struct SmallStrainMechanicalBehaviourIntegratorCRTPBase<
        SmallStrainMechanicalBehaviourIntegrator<H>>;

    /*!
     * \brief update the strain with the contribution of the given node
     * \param[in] g: strain
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateStrain(real *const,
                      const mfem::Vector &,
                      const mfem::DenseMatrix &,
                      const size_type);
    /*!
     * \brief update the inner forces of the given node  with the
     * contribution of the stress of an integration point.
     *
     * \param[out] Fe: inner forces
     * \param[in] s: stress
     * \param[in] dshape: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateInnerForces(mfem::Vector &,
                           const mgis::span<const real>&,
                           const mfem::DenseMatrix &,
                           const real,
                           const size_type) const;
    /*!
     * \brief update the stiffness matrix of the given node with the
     * contribution of the consistent tangent operator of an
     * integration point.
     *
     * \param[out] Ke: inner forces
     * \param[in] Kip: stress
     * \param[in] dshape: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &,
                               const mgis::span<const real>&,
                               const mfem::DenseMatrix &,
                               const real,
                               const size_type) const;

  };  // end of struct SmallStrainMechanicalBehaviourIntegrator

}  // end of namespace mfem_mgis

#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.ixx"

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX */
