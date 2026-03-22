#ifndef LIB_MFEM_MGIS_ISOTROPICPLANESTRESSBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ISOTROPICPLANESTRESSBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/FBarBehaviourIntegratorCRTPBase.hxx"
#include "MFEMMGIS/PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct FBarIsotropicPlaneStressBehaviourIntegrator;

  /*!
   * \brief partial specialisation of the `BehaviourIntegratorTraits`  * class
   * for the
   * `FBarIsotropicPlaneStressBehaviourIntegrator`
   * behaviour integrator */
  template <>
  struct BehaviourIntegratorTraits<
      FBarIsotropicPlaneStressBehaviourIntegrator> {
    static constexpr size_type unknownsSize = 2;
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   */
  struct MFEM_MGIS_EXPORT FBarIsotropicPlaneStressBehaviourIntegrator
      : FBarBehaviourIntegratorCRTPBase<
            FBarIsotropicPlaneStressBehaviourIntegrator,
            Hypothesis::PLANESTRESS>,
        PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase {
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
    FBarIsotropicPlaneStressBehaviourIntegrator(
        const FiniteElementDiscretization &,
        const size_type,
        std::unique_ptr<const Behaviour>);
    /*!
     * \return the rotation matrix associated with the given integration
     * point
     * \param[in] i: integration points
     */
    inline RotationMatrix getRotationMatrix(const size_type) const;

    inline void rotateGradients(std::span<real>, const RotationMatrix &);

    inline std::span<const real> rotateThermodynamicForces(
        std::span<const real>, const RotationMatrix &);

    inline void rotateTangentOperatorBlocks(std::span<real>,
                                            const RotationMatrix &);
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
    [[nodiscard]] bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept override;

    //! \brief destructor
    ~FBarIsotropicPlaneStressBehaviourIntegrator() override;

   protected:
    //! \brief allow the CRTP base class the protected members
    friend struct FBarBehaviourIntegratorCRTPBase<
        FBarIsotropicPlaneStressBehaviourIntegrator,
        Hypothesis::PLANESTRESS>;
    /*!
     * \return the integration rule for the given element and element
     * transformation.
     *
     * \param[in] e: element
     * \param[in] tr: element transformation
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
  };  // end of struct
      // FBarIsotropicPlaneStressBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ISOTROPICPLANESTRESSBEHAVIOURINTEGRATOR_HXX*/
