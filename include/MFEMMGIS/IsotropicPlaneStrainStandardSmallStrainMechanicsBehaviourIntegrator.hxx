#ifndef LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDSMALLSTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDSMALLSTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator;

  /*!
   * \brief specialisation of the `BehaviourIntegratorTraits` class for the
   * `IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator`
   * behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<
      IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator> {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 2;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   * \brief behaviour integrator for isotropic small strain mechanical
   * behaviours under the plane strain hypothesis
   */
  struct MFEM_MGIS_EXPORT
      IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator final
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator> {
    /*!
     * \brief a constant value used for the computation of
     * symmetric tensors
     */
    static constexpr const auto icste = real{0.70710678118654752440};
    //! \brief a dummy structure
    struct RotationMatrix {};
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator(
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
    ~IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator()
        override;

   protected:
    //! \brief allow the CRTP base class to access the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>;
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
     * \brief update the strain with the contribution of the
     * given node
     * \param[in, out] g: strain
     * \param[in] u: nodal displacements
     * \param[in] dN: derivatives of the shape functions
     * \param[in] ni: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &u,
                         const mfem::DenseMatrix &dN,
                         const size_type ni) noexcept;
    /*!
     * \brief update the inner forces of the given node with
     * the contribution of the stress of an integration point.
     *
     * \param[in, out] Fe: inner forces
     * \param[in] s: stress
     * \param[in] dN: derivatives of the shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateInnerForces(mfem::Vector &Fe,
                           const std::span<const real> &s,
                           const mfem::DenseMatrix &dN,
                           const real w,
                           const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of an
     * integration point.
     *
     * \param[in, out] Ke: stiffness matrix
     * \param[in] Kip: consistent tangent operator of the integration point
     * \param[in] dN: derivatives of the shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::DenseMatrix &dN,
                               const real w,
                               const size_type ni) const noexcept;

  };  // end of struct
      // IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDSMALLSTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX*/
