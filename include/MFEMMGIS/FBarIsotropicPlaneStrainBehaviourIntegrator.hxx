/*!
 * \file   include/MFEMMGIS/FBarIsotropicPlaneStrainBehaviourIntegrator.hxx
 * \brief  This file declares the `FBarIsotropicPlaneStrainBehaviourIntegrator`
 * class
 */

#ifndef LIB_MFEM_MGIS_ISOTROPICPLANESTRAINBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ISOTROPICPLANESTRAINBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/FBarBehaviourIntegratorCRTPBase.hxx"
#include "MFEMMGIS/PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct FBarIsotropicPlaneStrainBehaviourIntegrator;

  /*!
   * \brief specialisation of the `BehaviourIntegratorTraits` class for the
   * `FBarIsotropicPlaneStrainBehaviourIntegrator` behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<
      FBarIsotropicPlaneStrainBehaviourIntegrator> {
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
   * \brief FBar behaviour integrator for isotropic finite strain mechanical
   * behaviours under the plane strain hypothesis
   */
  struct MFEM_MGIS_EXPORT FBarIsotropicPlaneStrainBehaviourIntegrator
      : FBarBehaviourIntegratorCRTPBase<
            FBarIsotropicPlaneStrainBehaviourIntegrator,
            Hypothesis::PLANESTRAIN>,
        PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase {
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
    FBarIsotropicPlaneStrainBehaviourIntegrator(
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
    //
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
    /*!
     * \return true, the jacobian depends on the current solution through the
     * deformation gradient at the center of the element
     */
    [[nodiscard]] bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept override;

    //! \brief destructor
    ~FBarIsotropicPlaneStrainBehaviourIntegrator() override;

   protected:
    //! \brief allow the CRTP base class to access the protected members
    friend struct FBarBehaviourIntegratorCRTPBase<
        FBarIsotropicPlaneStrainBehaviourIntegrator,
        Hypothesis::PLANESTRAIN>;
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
  };  // end of struct
      // FBarIsotropicPlaneStrainBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ISOTROPICPLANESTRAINBEHAVIOURINTEGRATOR_HXX*/
