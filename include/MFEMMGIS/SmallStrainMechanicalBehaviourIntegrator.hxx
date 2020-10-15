/*!
 * \file   include/MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#ifndef LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX

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
    //! \brief destructor
    ~SmallStrainMechanicalBehaviourIntegratorBase() override;
  };  // end of struct SmallStrainMechanicalBehaviourIntegratorBase

  /*!
   * \brief behaviour integrators based on small strain
   * mechanical behaviours.
   *
   * \param[in] H: modelling hypothesis
   */
  template <Hypothesis H>
  struct MFEM_MGIS_EXPORT SmallStrainMechanicalBehaviourIntegrator
      : SmallStrainMechanicalBehaviourIntegratorBase {
    /*!
     * \brief constructor
     * \param[in] fs: finite element space.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    SmallStrainMechanicalBehaviourIntegrator(const mfem::FiniteElementSpace &,
                                             const size_type,
                                             std::unique_ptr<const Behaviour>);

    void computeInnerForces(const mfem::FiniteElement &,
                            mfem::ElementTransformation &,
                            const mfem::Vector &,
                            mfem::Vector &) override;
    void computeStiffnessMatrix(const mfem::FiniteElement &,
                                mfem::ElementTransformation &,
                                const mfem::Vector &,
                                mfem::DenseMatrix &) override;
    //! \brief destructor
    ~SmallStrainMechanicalBehaviourIntegrator() override;

  };  // end of struct SmallStrainMechanicalBehaviourIntegrator

}  // end of namespace mfem_mgis

#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.ixx"

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX */
