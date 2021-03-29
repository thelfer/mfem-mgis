/*!
 * \file   MFEMMGIS/BehaviourIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   27/08/2020
 */

#ifndef LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX
#define LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX

#include <memory>
#include "MGIS/Span.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Material;

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
  struct MFEM_MGIS_EXPORT BehaviourIntegrator {
    /*!
     * \brief set the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real) = 0;
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
     * \brief method call at the beginning of each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real) = 0;
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
    //! \return the underlying material
    virtual Material &getMaterial() = 0;
    //! \return the underlying material
    virtual const Material &getMaterial() const = 0;
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(mgis::span<const real>) = 0;
    //! \brief destructor
    virtual ~BehaviourIntegrator();
  };  // end of struct BehaviourIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BEHAVIOURINTEGRATOR_HXX */
