#ifndef LIB_MFEM_MGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX
#define LIB_MFEM_MGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX

#include <span>
#include <array>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  struct MFEM_MGIS_EXPORT
      PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase {
   protected:
    /*!
     * \brief update the deformation gradient with the contribution of the
     * given node
     * \param[in] g: deformation gradient
     * \param[in] u: nodal displacements
     * \param[in] dshape: derivatives of the shape function
     * \param[in] n: node index
     */
    void updateGradients(std::span<real> &,
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
                           const std::span<const real> &,
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
     * \param[in] dN: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &,
                               const std::span<const real> &,
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
     * \param[in] dN1: derivatives of the shape function
     * \param[in] dN2: derivatives of the shape function
     * \param[in] w: weight of the integration point
     * \param[in] n: node index
     */
    void updateStiffnessMatrix(mfem::DenseMatrix &,
                               const std::span<const real> &,
                               const mfem::DenseMatrix &,
                               const mfem::DenseMatrix &,
                               const real,
                               const size_type) const noexcept;
  };  // end of struct
      // PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase

}  // end of namespace mfem_mgis

#include "MFEMMGIS/PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase.ixx"

#endif /* LIB_MFEM_MGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_HXX*/
