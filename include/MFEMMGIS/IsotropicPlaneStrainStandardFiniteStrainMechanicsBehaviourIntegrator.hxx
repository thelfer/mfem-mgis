#ifndef LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"
#include "MFEMMGIS/PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator;

  /*!
   * \brief partial specialisation of the `BehaviourIntegratorTraits`  * class
   * for the
   * `IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator`
   * behaviour integrator */
  template <>
  struct BehaviourIntegratorTraits<
      IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator> {
    static constexpr size_type unknownsSize = 2;
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    static constexpr bool updateExternalStateVariablesFromUnknownsValues =
        false;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   */
  struct MFEM_MGIS_EXPORT
      IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>,
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
    IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator(
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
    ~IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator()
        override;

   protected:
    //! \brief allow the CRTP base class the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>;
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
  };  // end of struct
      // IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ISOTROPICPLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATOR_HXX*/
