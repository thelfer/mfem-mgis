/*!
 * \file   MFEMMGIS/AbstractBehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX

#include <array>
#include <memory>
#include <span>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Material;
  struct PartialQuadratureSpace;
  struct ImmutablePartialQuadratureFunctionView;
  enum struct IntegrationType;

  /*!
   * \brief abstract class for all behaviour integrators
   *
   * This class provides methods to:
   *
   * - integrate the behaviour over the time step
   * - compute the nodal forces du to the material reaction (see the
   *   `updateResidual` method).
   * - compute the stiffness matrix (see the `updateJacobian` method).
   */
  struct MFEM_MGIS_EXPORT AbstractBehaviourIntegrator {
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real) = 0;
    //! \return the partial quadrature space
    virtual const PartialQuadratureSpace &getPartialQuadratureSpace() const = 0;
    /*!
     * \return the integration rule for the given element and element
     * transformation
     * \param[in] e: element
     * \param[in] tr: element transformation
     */
    virtual const mfem::IntegrationRule &getIntegrationRule(
        const mfem::FiniteElement &,
        const mfem::ElementTransformation &) const = 0;
    /*!
     * \brief return the weight of the integration point, taking the
     * modelling hypothesis into account
     * \param[in] tr: element transformation
     * \param[in] ip: integration point
     */
    virtual real getIntegrationPointWeight(
        mfem::ElementTransformation &,
        const mfem::IntegrationPoint &) const = 0;
    /*!
     * \brief method call at the beginning of each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real) = 0;
    /*!
     * \brief integrate the mechanical behaviour over the time step
     * If successful, the value of the stress, consistent tangent
     * operator and internal state variables are updated.
     *
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     */
    virtual bool integrate(const mfem::FiniteElement &,
                           mfem::ElementTransformation &,
                           const mfem::Vector &,
                           const IntegrationType) = 0;
    /*!
     * \brief compute the contribution of the given element to the inner forces
     * \param[out] Fe: inner forces
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     */
    virtual void computeInnerForces(mfem::Vector &,
                                    const mfem::FiniteElement &,
                                    mfem::ElementTransformation &) = 0;
    /*!
     * \brief compute the contribution of the given element to the residual
     * \param[out] Fe: element contribution to the residual
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     */
    virtual void updateResidual(mfem::Vector &,
                                const mfem::FiniteElement &,
                                mfem::ElementTransformation &,
                                const mfem::Vector &) = 0;
    /*!
     * \brief compute the contribution of the given element to the jacobian
     * \param[out] Ke: element stiffness matrix
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimation of the displacement field
     */
    virtual void updateJacobian(mfem::DenseMatrix &,
                                const mfem::FiniteElement &,
                                mfem::ElementTransformation &,
                                const mfem::Vector &) = 0;
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
    virtual bool hasMaterial() const = 0;
    //! \return the underlying material
    virtual Material &getMaterial() = 0;
    //! \return the underlying material
    virtual const Material &getMaterial() const = 0;
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(std::span<const real>) = 0;
    /*!
     * \return if the current solution is required for assembling the residual
     * \note this is required when creating a linear operator evaluating the
     * the residual
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForResidualAssembly()
        const noexcept = 0;
    /*!
     * \return if the current solution is required for assembling the jacobian
     * \note this is required when creating a linear operator evaluating the
     * the residual
     */
    [[nodiscard]] virtual bool requiresCurrentSolutionForJacobianAssembly()
        const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractBehaviourIntegrator();
  };  // end of struct AbstractBehaviourIntegrator

  /*!
   * \return the measure of the mesh (volume in 3D and axisymmetrical
   * hypotheses, surface in other bidimensional hypotheses) on which is built
   * the behaviour integrator. \param[in] bi: behaviour integrator
   */
  MFEM_MGIS_EXPORT real computeMeasure(const AbstractBehaviourIntegrator &);

  /*!
   * \return the integral of a partial quadrature function
   * \param[in] f: function
   */
  template <typename ValueType>
  ValueType computeIntegral(const AbstractBehaviourIntegrator &,
                            const ImmutablePartialQuadratureFunctionView &);
  /*!
   * \return the integral of a partial quadrature function
   * \param[in] bi: behaviour integrator
   * \param[in] f: function
   */
  template <>
  MFEM_MGIS_EXPORT real
  computeIntegral(const AbstractBehaviourIntegrator &,
                  const ImmutablePartialQuadratureFunctionView &);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTBEHAVIOURINTEGRATOR_HXX */
