/*!
 * \file   MFEMMGIS/AbstractBehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX

#include <span>
#include <array>
#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Material;
  struct PartialQuadratureSpace;
  struct ImmutablePartialQuadratureFunctionView;
  struct AbstractQPEvaluator;
  enum struct IntegrationType;

  /*!
   * \brief abstract class for all behaviour integrators
   *
   * This class provides methods to:
   *
   * - integrate the behaviour over the time step
   * - compute the nodal forces due to the material reaction (see the
   *   `updateResidual` method).
   * - compute the stiffness matrix (see the `updateJacobian` method).
   */
  struct MFEM_MGIS_EXPORT AbstractBehaviourIntegrator {
    //! \return the current time increment
    virtual real getTimeIncrement() const noexcept = 0;
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real dt) = 0;
    //! \return the partial quadrature space
    virtual const PartialQuadratureSpace &getPartialQuadratureSpace()
        const noexcept = 0;
    /*!
     * \brief return the integration rule for the given element and element
     * transformation
     * \param[in] e: element
     * \param[in] tr: element transformation
     * \return the integration rule
     */
    virtual const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &e,
        const mfem::ElementTransformation &tr) const = 0;
    /*!
     * \brief return the weight of the integration point, taking the
     * modelling hypothesis into account
     * \param[in] tr: element transformation
     * \param[in] ip: integration point
     * \return the weight of the integration point
     */
    virtual real getIntegrationPointWeight(
        mfem::ElementTransformation &tr,
        const mfem::IntegrationPoint &ip) const = 0;
    /*!
     * \brief method called at the beginning of each resolution
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool setup(Context &ctx,
                                     const real t,
                                     const real dt) noexcept = 0;
    /*!
     * \brief integrate the behaviour over the time step
     *
     * If successful, the values of the thermodynamic forces, consistent
     * tangent operator and internal state variables are updated.
     *
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     * \return true on success
     */
    [[nodiscard]] virtual bool integrate(const mfem::FiniteElement &e,
                                         mfem::ElementTransformation &tr,
                                         const mfem::Vector &u,
                                         const IntegrationType it) = 0;
    /*!
     * \brief compute the contribution of the given element to the inner forces
     * \param[out] Fe: inner forces
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     */
    virtual void computeInnerForces(mfem::Vector &Fe,
                                    const mfem::FiniteElement &e,
                                    mfem::ElementTransformation &tr) = 0;
    /*!
     * \brief compute the contribution of the given element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    virtual void updateResidual(mfem::Vector &Fe,
                                const mfem::FiniteElement &e,
                                mfem::ElementTransformation &tr,
                                const mfem::Vector &u) = 0;
    /*!
     * \brief compute the contribution of the given element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in, out] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     */
    virtual void updateJacobian(mfem::DenseMatrix &Ke,
                                const mfem::FiniteElement &e,
                                mfem::ElementTransformation &tr,
                                const mfem::Vector &u) = 0;
    /*!
     * \brief clean-up at the end of a resolution
     *
     * \param[in, out] ctx: execution context
     * \return true on success
     *
     * \note this method is the counterpart of the `setup` method and is meant
     * to release any memory allocated by this method.
     */
    [[nodiscard]] virtual bool cleanup(Context &ctx) noexcept = 0;
    /*!
     * \brief revert the internal state variables.
     *
     * The values of the internal state variables at the beginning of the time
     * step are copied on the values of the internal state variables at
     * end of the time step.
     */
    virtual void revert() = 0;
    /*!
     * \brief update the internal state variables.
     *
     * The values of the internal state variables at the end of the time step
     * are copied on the values of the internal state variables at beginning of
     * the time step.
     */
    virtual void update() = 0;
    /*!
     * \return if the call to getMaterial is valid
     *
     * This has been introduced to be able to build behaviour integrators not
     * built on MGIS and MFront.
     */
    virtual bool hasMaterial() const noexcept = 0;
    /*!
     * \brief return the underlying material
     * \param[in, out] ctx: execution context
     * \return the underlying material, if any
     */
    [[nodiscard]] virtual OptionalReference<Material> getMaterial(
        Context &ctx) noexcept = 0;
    /*!
     * \brief return the underlying material
     * \param[in, out] ctx: execution context
     * \return the underlying material, if any
     */
    [[nodiscard]] virtual OptionalReference<const Material> getMaterial(
        Context &ctx) const noexcept = 0;
    /*!
     * \brief set the macroscopic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(std::span<const real> g) = 0;
    /*!
     * \return if the current solution is required for assembling the residual
     * \note this is required when creating a linear operator evaluating the
     * residual
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForResidualAssembly()
        const noexcept = 0;
    /*!
     * \return if the current solution is required for assembling the jacobian
     * \note this is required when creating a linear operator evaluating the
     * jacobian
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept = 0;
    /*!
     * \brief set the value of a material property
     *
     * \param[in, out] ctx: execution context
     * \param[in] name: name of the material property
     * \param[in] e: evaluator of the material property
     * \param[in] ts: time step stage
     * \return true on success
     */
    [[nodiscard]] virtual bool setMaterialProperty(
        Context &ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept = 0;
    /*!
     * \brief set the value of an external state variable
     *
     * \param[in, out] ctx: execution context
     * \param[in] name: name of the external state variable
     * \param[in] e: evaluator of the external state variable
     * \param[in] ts: time step stage
     * \return true on success
     */
    [[nodiscard]] virtual bool setExternalStateVariable(
        Context &ctx,
        std::string_view name,
        std::shared_ptr<const AbstractQPEvaluator> e,
        const TimeStepStage ts) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractBehaviourIntegrator();
  };  // end of struct AbstractBehaviourIntegrator

  /*!
   * \brief return the measure of the mesh on which the behaviour integrator is
   * built: volume in 3D and axisymmetrical hypotheses, surface in other
   * bidimensional hypotheses.
   * \param[in] bi: behaviour integrator
   * \return the measure of the mesh
   */
  MFEM_MGIS_EXPORT real computeMeasure(const AbstractBehaviourIntegrator &bi);

  /*!
   * \brief return the integral of a partial quadrature function
   * \param[in] bi: behaviour integrator
   * \param[in] f: function
   * \return the integral of the function
   */
  template <typename ValueType>
  ValueType computeIntegral(const AbstractBehaviourIntegrator &bi,
                            const ImmutablePartialQuadratureFunctionView &f);
  /*!
   * \return the integral of a partial quadrature function
   * \param[in] bi: behaviour integrator
   * \param[in] f: function
   */
  template <>
  MFEM_MGIS_EXPORT real
  computeIntegral(const AbstractBehaviourIntegrator &bi,
                  const ImmutablePartialQuadratureFunctionView &f);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX */
