/*!
 * \file
 * \brief  This file declares the
 * `OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator`
 * class
 */

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
   * \brief specialisation of the `BehaviourIntegratorTraits` class for the
   * `OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator`
   * behaviour integrator
   */
  template <>
  struct BehaviourIntegratorTraits<
      OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator> {
    //! \brief number of components of the unknowns
    static constexpr size_type unknownsSize = 1;
    //! \brief if the computation of the gradients requires the shape functions
    static constexpr bool gradientsComputationRequiresShapeFunctions = false;
    /*!
     * \brief if the computation of the gradients requires the derivatives of
     * the shape functions
     */
    static constexpr bool
        gradientsComputationRequiresShapeFunctionsDerivatives = true;
    //! \brief if the external state variables are updated from the unknowns
    static constexpr bool updateExternalStateVariablesFromUnknownsValues = true;
  };  // end of struct BehaviourIntegratorTraits<>

  /*!
   * \brief behaviour integrator for orthotropic stationary non linear heat
   * transfer behaviours under the plane strain hypothesis
   */
  struct MFEM_MGIS_EXPORT
      OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator
          final
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator> {
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
    /*!
     * \brief method called at the beginning of each resolution
     *
     * In addition to the base class setup, store a pointer to the values of
     * the `Temperature` external state variable, which must be defined,
     * mutable and not uniform.
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
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
     */
    inline void rotateGradients(std::span<real> g, const RotationMatrix &r);

    /*!
     * \brief rotate the thermodynamic forces in the global frame
     * \param[in] s: thermodynamic forces in the material frame
     * \param[in] r: rotation matrix
     * \return the thermodynamic forces in the global frame
     */
    inline std::array<real, 2> rotateThermodynamicForces(
        std::span<const real> s, const RotationMatrix &r);

    /*!
     * \brief rotate the tangent operator blocks in the global frame
     * \param[in, out] Kip: tangent operator blocks
     * \param[in] r: rotation matrix
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
    ~OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator()
        override;

   protected:
    //! \brief allow the CRTP base class to access the protected members
    friend struct StandardBehaviourIntegratorCRTPBase<
        OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator>;
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
     * \brief update the temperature gradient with the contribution of the
     * given node
     * \param[in, out] g: temperature gradient
     * \param[in] u: nodal temperatures
     * \param[in] dN: derivatives of the shape functions
     * \param[in] ni: node index
     */
    void updateGradients(std::span<real> &g,
                         const mfem::Vector &u,
                         const mfem::DenseMatrix &dN,
                         const size_type ni) noexcept;
    /*!
     * \brief update the inner forces of the given node with
     * the contribution of the heat flux of an integration point.
     *
     * \param[in, out] Fe: inner forces
     * \param[in] s: heat flux
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
     * \param[in] N: shape function values
     * \param[in] dN: derivatives of the shape functions
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                               const std::span<const real> &Kip,
                               const mfem::Vector &N,
                               const mfem::DenseMatrix &dN,
                               const real w,
                               const size_type ni) const noexcept;

    /*!
     * \brief pointer to the external state variable
     * associated with the unknown
     */
    real *uesv = nullptr;
  };  // end of struct
      // OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ORTHOTROPICPLANESTRAINSTATIONARYNONLINEARHEATTRANSFERBEHAVIOURINTEGRATOR_HXX*/
