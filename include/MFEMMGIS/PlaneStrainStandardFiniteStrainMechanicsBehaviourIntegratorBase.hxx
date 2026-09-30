#ifndef LIB_MFEM_MGIS_PLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX
#define LIB_MFEM_MGIS_PLANESTRAINSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX

#include <span>
#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  struct MFEM_MGIS_EXPORT
      PlaneStrainStandardFiniteStrainMechanicsBehaviourIntegratorBase {
   protected:
    /*!
     * \brief update the deformation gradient with the contribution of the
     * given node
     * \param[in] g: deformation gradient
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateGradients(std::span<real>& g,
                         const mfem::Vector& u,
                         const mfem::DenseMatrix& dN,
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
    void updateInnerForces(mfem::Vector& Fe,
                           const std::span<const real>& s,
                           const mfem::DenseMatrix& dN,
                           const real w,
                           const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of  * an
     * integration point.
     *
     * \param[out] Ke: inner forces
     * \param[in] Kip: stress
     * \param[in] dN: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix& Ke,
                               const std::span<const real>& Kip,
                               const mfem::DenseMatrix& dN,
                               const real w,
                               const size_type ni) const noexcept;
    /*!
     * \brief update the stiffness matrix of the given node
     * with the contribution of the consistent tangent operator of  * an
     * integration point.
     *
     * \param[out] Ke: inner forces
     * \param[in] Kip: stress
     * \param[in] dN1: derivatives of the shape function
     * \param[in] dN2: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
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
