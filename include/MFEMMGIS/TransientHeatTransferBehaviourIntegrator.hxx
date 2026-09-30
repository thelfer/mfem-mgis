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
   * \brief partial specialisation of the `BehaviourIntegratorTraits`
   * class for the `TransientHeatTransferBehaviourIntegrator`
   * behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<TransientHeatTransferBehaviourIntegrator> {
    static constexpr size_type unknownsSize = 1;
    static constexpr bool gradientsComputationRequiresShapeFunctions = true;
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = false;
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
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
     * \return the rotation matrix associated with the given integration
     * point
     * \param[in] i: integration points
     */
    inline RotationMatrix getRotationMatrix(const size_type i) const;

    inline void rotateGradients(std::span<real> g, const RotationMatrix &r);

    inline std::span<const real> rotateThermodynamicForces(
        std::span<const real> s, const RotationMatrix &r);

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
    //! \brief allow the CRTP base class the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        TransientHeatTransferBehaviourIntegrator>;
    /*!
     * \return the integration rule for the given element and  * element
     * transformation. \param[in] e: element \param[in] tr: element
     * transformation
     */
    static const mfem::IntegrationRule &selectIntegrationRule(
        const mfem::FiniteElement &e, const mfem::ElementTransformation &t);
    /*!
     * \brief build the quadrature space for the given  * material
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     */
    static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
        const FiniteElementDiscretization &fed, const size_type m);
    /*!
     * \brief update the strain with the contribution of the
     * given node
     * \param[in] g: strain
     * \param[in] u: nodal displacements
     * \param[in] N: shape functions
     * \param[in] n: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &T,
                         const mfem::Vector &N,
                         const size_type ni) noexcept;
    /*!
     * \brief update the inner forces of the given node  with
     * the contribution of the stress of an integration point.
     *
     * \param[out] Fe: inner forces
     * \param[in] s: stress
     * \param[in] N: shape functions
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateInnerForces(mfem::Vector &Fe,
                           const std::span<const real> &s,
                           const mfem::Vector &N,
                           const real w,
                           const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of  * an
     * integration point.
     *
     * \param[out] Ke: inner forces
     * \param[in] Kip: stress
     * \param[in] N: shape functions
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::Vector &N,
                               const real w,
                               const size_type ni) const noexcept;
  };  // end of struct TransientHeatTransferBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_TRANSIENT_HEAT_TRANSFERT_BEHAVIOUR_INTEGRATOR_HXX*/
