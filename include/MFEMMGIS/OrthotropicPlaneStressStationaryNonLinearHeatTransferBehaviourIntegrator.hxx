#ifndef LIB_MFEM_MGIS_ORTHOTROPICPLANESTRESSSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ORTHOTROPICPLANESTRESSSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"
#include "MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  // forward declaration
  struct
      OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator;

  /*!
   * \brief partial specialisation of the `BehaviourIntegratorTraits`  * class
   * for the
   * `OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator`
   * behaviour integrator */
  template <>
  struct BehaviourIntegratorTraits<
      OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator> {
    //! \brief size of the unknowns
    static constexpr size_type unknownsSize = 1;
    //! \brief
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    //! \brief
    static constexpr bool updateExternalStateVariablesFromUnknownsValues = true;
  };  // end of struct
      // BehaviourIntegratorTraits<OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>

  /*!
   */
  struct MFEM_MGIS_EXPORT
      OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator
          final
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator> {
    /*!
     * \brief a constant value used for the computation of
     * symmetric tensors
     */
    static constexpr const auto icste = real{0.70710678118654752440};

    //! \brief a simple alias
    using RotationMatrix =
        std::array<real, 9u>; /*!
                               * \brief constructor
                               * \param[in] fed: finite element discretization.
                               * \param[in] m: material attribute.
                               * \param[in] b_ptr: behaviour
                               */
    OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator(
        const FiniteElementDiscretization &,
        const size_type,
        std::unique_ptr<const Behaviour>);

    void setup(const real, const real) override;
    /*!
     * \brief update the external state variable
     * corresponding to the unknowns.
     *
     * \param[in] u: unknowns
     * \param[in] N: values of the shape functions
     * \param[in] o: offset associated with the
     * current integration point
     */
    void updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &,
                                                        const mfem::Vector &,
                                                        const size_type);

    /*!
     * \return the rotation matrix associated with the given  * integration
     * point \param[in] i: integration points
     */
    inline RotationMatrix getRotationMatrix(const size_type) const;

    inline void rotateGradients(mgis::span<real>, const RotationMatrix &);

    inline std::array<real, 2> rotateThermodynamicForces(
        mgis::span<const real>, const RotationMatrix &);

    inline void rotateTangentOperatorBlocks(mgis::span<real>,
                                            const RotationMatrix &);

    const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &,
        const mfem::ElementTransformation &) const override;

    real getIntegrationPointWeight(mfem::ElementTransformation &,
                                   const mfem::IntegrationPoint &) const
        noexcept override;

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

    //! \brief destructor
    ~OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator()
        override;

   protected:
    //! \brief allow the CRTP base class the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>;
    /*!
     * \return the integration rule for the given element and  * element
     * transformation. \param[in] e: element \param[in] tr: element
     * transformation
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
    /*!
     * \brief update the strain with the contribution of the
     * given node
     * \param[in] g: strain
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateGradients(mgis::span<real> &,
                         const mfem::Vector &,
                         const mfem::DenseMatrix &,
                         const size_type) noexcept;
    /*!
     * \brief update the inner forces of the given node  with
     * the contribution of the stress of an integration point.
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
                           const size_type) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of  * an
     * integration point.
     *
     * \param[out] Ke: inner forces
     * \param[in] Kip: stress
     * \param[in] N: shape function values
     * \param[in] dN: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &,
                               const mgis::span<const real> &,
                               const mfem::Vector &,
                               const mfem::DenseMatrix &,
                               const real,
                               const size_type) const noexcept;

    //! rief the rotation matrix
    RotationMatrix2D rotation_matrix;

    /*!
     * \brief pointer to the external state variable
     * associated with the unknown
     */
    real *uesv = nullptr;
  };  // end of struct
      // OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ORTHOTROPICPLANESTRESSSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX*/
