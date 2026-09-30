/*!
 * \file
 * \brief  This file declares the
 * `PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase` class
 */

#ifndef LIB_MFEM_MGIS_PLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX
#define LIB_MFEM_MGIS_PLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX

#include <span>
#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class of the behaviour integrators for finite strain
   * mechanical behaviours under the plane strain hypothesis
   */
  struct MFEM_MGIS_EXPORT
      PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase {
   protected:
    /*!
     * \brief update the deformation gradient with the contribution of the
     * given node
     * \param[in, out] g: deformation gradient
     * \param[in] u: nodal displacements
     * \param[in] dN: derivatives of the shape functions
     * \param[in] ni: node index
     */
    void updateGradients(std::span<real>& g,
                         const mfem::Vector& u,
                         const mfem::DenseMatrix& dN,
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
    void updateInnerForces(mfem::Vector& Fe,
                           const std::span<const real>& s,
                           const mfem::DenseMatrix& dN,
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
    void updateStiffnessMatrix(mfem::DenseMatrix& Ke,
                               const std::span<const real>& Kip,
                               const mfem::DenseMatrix& dN,
                               const real w,
                               const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of an
     * integration point.
     *
     * \param[in, out] Ke: stiffness matrix
     * \param[in] Kip: consistent tangent operator of the integration point
     * \param[in] dN1: derivatives of the shape functions, associated with the
     * rows of the stiffness matrix
     * \param[in] dN2: derivatives of the shape functions, associated with the
     * columns of the stiffness matrix
     * \param[in] w: weight of the integration point
     * \param[in] ni: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix& Ke,
                               const std::span<const real>& Kip,
                               const mfem::DenseMatrix& dN1,
                               const mfem::DenseMatrix& dN2,
                               const real w,
                               const size_type ni) const noexcept;
  };

}  // end of namespace mfem_mgis

#include "MFEMMGIS/PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase.ixx"

#endif /* LIB_MFEM_MGIS_PLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX*/
