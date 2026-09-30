#ifndef LIB_MFEM_MGIS_ORTHOTROPICPLANESTRAINSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ORTHOTROPICPLANESTRAINSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX

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
      OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator;

  /*!
   * \brief partial specialisation of the `BehaviourIntegratorTraits`  * class
   * for the
   * `OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator`
   * behaviour integrator */
  template <>
  struct BehaviourIntegratorTraits<
      OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator> {
    static constexpr size_type unknownsSize = 1;
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    static constexpr bool updateExternalStateVariablesFromUnknownsValues = true;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   */
  struct MFEM_MGIS_EXPORT
      OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator
          final
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator> {
    /*!
     * \brief a constant value used for the computation of
     * symmetric tensors
     */
    static constexpr const auto icste = real{0.70710678118654752440};
    //! \brief a simple alias
    using RotationMatrix = std::array<real, 9u>;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b_ptr: behaviour
     */
    OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator(
        const FiniteElementDiscretization &fed,
        const size_type m,
        std::unique_ptr<const Behaviour> b_ptr);
    //
    [[nodiscard]] bool setup(Context &ctx,
                             const real t,
                             const real dt) noexcept override;
    /*!
     * \brief update the external state variable
     * corresponding to the unknowns.
     *
     * \param[in] u: unknowns
     * \param[in] N: values of the shape functions
     * \param[in] o: offset associated with the
     * current integration point
     */
    void updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &u,
                                                        const mfem::Vector &N,
                                                        const size_type o);

    /*!
     * \return the rotation matrix associated with the given integration
     * point
     * \param[in] i: integration points
     */
    inline RotationMatrix getRotationMatrix(const size_type i) const;

    inline void rotateGradients(std::span<real> g, const RotationMatrix &r);

    inline std::array<real, 2> rotateThermodynamicForces(
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
    ~OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator()
        override;

   protected:
    //! \brief allow the CRTP base class the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator>;
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
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &u,
                         const mfem::DenseMatrix &dN,
                         const size_type ni) noexcept;
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
    void updateInnerForces(mfem::Vector &Fe,
                           const std::span<const real> &s,
                           const mfem::DenseMatrix &dN,
                           const real w,
                           const size_type ni) const noexcept;
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
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::Vector &N,
                               const mfem::DenseMatrix &dN,
                               const real w,
                               const size_type ni) const noexcept;

    //! \brief the rotation matrix
    RotationMatrix2D rotation_matrix;

    /*!
     * \brief pointer to the external state variable
     * associated with the unknown
     */
    real *uesv = nullptr;
  };  // end of struct
      // OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ORTHOTROPICPLANESTRAINSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX*/
