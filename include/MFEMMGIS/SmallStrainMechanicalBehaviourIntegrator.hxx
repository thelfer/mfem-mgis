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
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief behaviour integrators based on small strain
   * mechanical behaviours.
   *
   * \param[in] H: modelling hypothesis
   */
  template <Hypothesis H>
  struct MFEM_MGIS_EXPORT SmallStrainMechanicalBehaviourIntegrator
      : StandardBehaviourIntegratorCRTPBase<
            SmallStrainMechanicalBehaviourIntegrator<H>> {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    SmallStrainMechanicalBehaviourIntegrator(const FiniteElementDiscretization &,
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
    friend struct StandardBehaviourIntegratorCRTPBase<
        SmallStrainMechanicalBehaviourIntegrator<H>>;
    /*!
     * \return the integration rule for the given element and element
     * transformation.
     * \param[in] e: element
     * \param[in] tr: element transformation
     */
    static const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &, const mfem::ElementTransformation &);
    /*!
     * \brief build the quadrature space for the given material
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        const FiniteElementDiscretization &, const size_type);
    /*!
     * \brief update the strain with the contribution of the given node
     * \param[in] g: strain
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateGradients(mgis::span<real> &,
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
                           const mgis::span<const real> &,
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
                               const mgis::span<const real> &,
                               const mfem::DenseMatrix &,
                               const real,
                               const size_type) const;

    };  // end of struct SmallStrainMechanicalBehaviourIntegrator

}  // end of namespace mfem_mgis

#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.ixx"

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_HXX */
