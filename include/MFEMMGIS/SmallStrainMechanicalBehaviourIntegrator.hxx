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
        std::unique_ptr<const Behaviour>);
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

   protected:
    /*!
     * \brief update the strain with the contribution of the given node
     * \param[in] g: strain
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    virtual void updateStrain(real *const,
                              const mfem::Vector &,
                              const mfem::DenseMatrix &,
                              const size_type) = 0;
    /*!
     * \brief update the inner forces the given node
     * \param[out] Fe: inner forces
     * \param[in] s: stress
     * \param[in] dshape: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    virtual void updateInnerForces(mfem::Vector &,
                                   const real *const,
                                   const mfem::DenseMatrix &,
                                   const real,
                                   const size_type) const = 0;

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
     * \brief This method is meant to implement the `computeInnerForces` method.
     */
    void implementComputeInnerForces(mfem::Vector &,
                                     const mfem::FiniteElement &,
                                     mfem::ElementTransformation &,
                                     const mfem::Vector &);
    /*!
     * \brief This method is meant to implement the `computeStiffnessMatrix`
     * method.
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
                                             std::unique_ptr<const Behaviour>);

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

    void updateStrain(real *const,
                      const mfem::Vector &,
                      const mfem::DenseMatrix &,
                      const size_type) override;
    void updateInnerForces(mfem::Vector &,
                           const real *const,
                           const mfem::DenseMatrix &,
                           const real,
                           const size_type) const override;

  };  // end of struct SmallStrainMechanicalBehaviourIntegrator

}  // end of namespace mfem_mgis

#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.ixx"

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX */
