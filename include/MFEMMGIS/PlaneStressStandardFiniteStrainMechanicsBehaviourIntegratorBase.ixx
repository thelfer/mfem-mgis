/*!
 * \file
 * MFEMMGIS/PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   19/03/2026
 */

#ifndef LIB_MFEMMGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_IXX
#define LIB_MFEMMGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_IXX

namespace mfem_mgis {

  inline void
  PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase::
      updateGradients(std::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto nnodes = dN.NumRows();
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    g[0] += dNi_0 * u_0;
    g[1] += dNi_1 * u_1;
    g[2] += 0;
    g[3] += dNi_1 * u_0;
    g[4] += dNi_0 * u_1;
  }  // end of updateGradients

  inline void
  PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase::
      updateInnerForces(mfem::Vector &Fe,
                        const std::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto nnodes = dN.NumRows();
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    Fe[ni_0] += w * (dNi_0 * s[0] + dNi_1 * s[3]);
    Fe[ni_1] += w * (s[4] * dNi_0 + s[1] * dNi_1);
  }  // end of updateInnerForces

  inline void
  PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase::
        updateStiffnessMatrix(Ke, Kip, dN, dN, w, ni);
  }  // end of updateStiffnessMatrix

  inline void
  PlaneStressStandardFiniteStrainMechanicsBehaviourIntegratorBase::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::DenseMatrix &dN1,
                            const mfem::DenseMatrix &dN2,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN1.NumRows();
    const auto dNi_0 = dN1(ni, 0);
    const auto dNi_1 = dN1(ni, 1);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN2(nj, 0);
      const auto dNj_1 = dN2(nj, 1);
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      Ke(ni_0, nj_0) += w * (dNi_1 * Kip[15] * dNj_0 + dNj_1 * dNi_0 * Kip[3] +
                             dNj_1 * dNi_1 * Kip[18] + Kip[0] * dNi_0 * dNj_0);
      Ke(ni_0, nj_1) += w * (Kip[16] * dNj_1 * dNi_1 + dNj_1 * Kip[1] * dNi_0 +
                             Kip[4] * dNi_0 * dNj_0 + Kip[19] * dNi_1 * dNj_0);
      Ke(ni_1, nj_0) += w * (dNj_1 * Kip[23] * dNi_0 + Kip[5] * dNi_1 * dNj_0 +
                             dNi_0 * Kip[20] * dNj_0 + dNj_1 * Kip[8] * dNi_1);
      Ke(ni_1, nj_1) += w * (dNi_1 * Kip[9] * dNj_0 + dNj_1 * dNi_0 * Kip[21] +
                             dNj_1 * dNi_1 * Kip[6] + dNi_0 * Kip[24] * dNj_0);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PLANESTRESSSTANDARDFINITESTRAINMECHANICSBEHAVIOURINTEGRATORBASE_IXX \
        */
