/*!
 * \file   include/MFEMMGIS/TransientHeatTransferBehaviourIntegrator.hxx
 * \brief  This file declares the `TransientHeatTransferBehaviourIntegrator`
 * class
 */

#ifndef LIB_MFEM_MGIS_TRANSIENT_HEAT_TRANSFERT_BEHAVIOUR_INTEGRATOR_HXX
#define LIB_MFEM_MGIS_TRANSIENT_HEAT_TRANSFERT_BEHAVIOUR_INTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct TransientHeatTransferBehaviourIntegrator;

  /*!
   * \brief specialisation of the `BehaviourIntegratorTraits`
   * class for the `TransientHeatTransferBehaviourIntegrator`
   * behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<TransientHeatTransferBehaviourIntegrator> {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 1;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = true;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = false;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   * \brief behaviour integrator for the transient term of the heat equation
   *
   * This behaviour integrator is meant to compute the inner forces associated
   * with the term:
   *
   * \f[
   * \frac{\partial h}{\partial T}\,\dot{T}
   * \f]
   *
   * where \f$h\f$ is the enthalpy per unit of volume.
   *
   * \note Most of the time, the derivative
   * \f$\frac{\partial h}{\partial T}\f$
   * is computed by the product of the mass density by the specific heat.
   */
  struct MFEM_MGIS_EXPORT TransientHeatTransferBehaviourIntegrator final
      : StandardBehaviourIntegratorCRTPBase<
            TransientHeatTransferBehaviourIntegrator> {
    //! \brief a dummy structure
    struct RotationMatrix {};
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    TransientHeatTransferBehaviourIntegrator(
        const FiniteElementDiscretization &fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr);
    /*!
     * \brief get the rotation matrix at the given integration point
     * \return the rotation matrix associated with the given integration
     * point
     * \param[in] i: offset of the integration point
     */
    inline RotationMatrix getRotationMatrix(const size_type i) const;

    /*!
     * \brief rotate the gradients in the material frame
     * \param[in, out] g: gradients
     * \param[in] r: rotation matrix
     * \note does nothing
     */
    inline void rotateGradients(std::span<real> g, const RotationMatrix &r);

    /*!
     * \brief rotate the thermodynamic forces in the global frame
     * \param[in] s: thermodynamic forces
     * \param[in] r: rotation matrix
     * \return the thermodynamic forces in the global frame
     * \note returns s unchanged
     */
    inline std::span<const real> rotateThermodynamicForces(
        std::span<const real> s, const RotationMatrix &r);

    /*!
     * \brief rotate the tangent operator blocks in the global frame
     * \param[in, out] Kip: tangent operator blocks
     * \param[in] r: rotation matrix
     * \note does nothing
     */
    inline void rotateTangentOperatorBlocks(std::span<real> Kip,
                                            const RotationMatrix &r);

    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const override;

    real getIntegrationPointWeight(
        mfem::ElementTransformation &tr,
        const mfem::IntegrationPoint &ip) const noexcept override;

    bool integrate(const mfem::FiniteElement &e,
                   mfem::ElementTransformation &tr,
                   const mfem::Vector &u,
                   const IntegrationType it) override;

    void updateResidual(mfem::Vector &Fe,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;

    void updateJacobian(mfem::DenseMatrix &Ke,
                        const mfem::FiniteElement &e,
                        mfem::ElementTransformation &tr,
                        const mfem::Vector &u) override;

    void computeInnerForces(mfem::Vector &Fe,
                            const mfem::FiniteElement &e,
                            mfem::ElementTransformation &tr) override;

    //! \brief destructor
    ~TransientHeatTransferBehaviourIntegrator() override;

   protected:
    //! \brief allow the CRTP base class to access the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        TransientHeatTransferBehaviourIntegrator>;
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
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \return the partial quadrature space
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        const FiniteElementDiscretization &fed, const size_type m);
    /*!
     * \brief update the temperature at the integration point with the
     * contribution of the given node
     * \param[in, out] g: temperature
     * \param[in] T: nodal temperatures
     * \param[in] N: shape functions
     * \param[in] ni: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &T,
                         const mfem::Vector &N,
                         const size_type ni) noexcept;
    /*!
     * \brief update the inner forces of the given node with
     * the contribution of the thermodynamic force of an integration point.
     *
     * \param[in, out] Fe: inner forces
     * \param[in] s: thermodynamic force
     * \param[in] N: shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateInnerForces(mfem::Vector &Fe,
                           const std::span<const real> &s,
                           const mfem::Vector &N,
                           const real w,
                           const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of an
     * integration point.
     *
     * \param[in, out] Ke: stiffness matrix
     * \param[in] Kip: consistent tangent operator of the integration point
     * \param[in] N: shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::Vector &N,
                               const real w,
                               const size_type ni) const noexcept;
  };  // end of struct TransientHeatTransferBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_TRANSIENT_HEAT_TRANSFERT_BEHAVIOUR_INTEGRATOR_HXX*/
