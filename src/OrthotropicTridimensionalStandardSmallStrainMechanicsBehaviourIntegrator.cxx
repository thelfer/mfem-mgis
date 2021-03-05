#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  inline void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_2_2 = dN(ni, 2);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto Bi_4_0 = dN(ni, 2) * icste;
    const auto Bi_4_2 = dN(ni, 0) * icste;
    const auto Bi_5_1 = dN(ni, 2) * icste;
    const auto Bi_5_2 = dN(ni, 1) * icste;
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    const auto u_2 = u[ni + 2 * nnodes];
    g[0] += u_0 * Bi_0_0;
    g[1] += u_1 * Bi_1_1;
    g[2] += u_2 * Bi_2_2;
    g[3] += u_0 * Bi_3_0 + u_1 * Bi_3_1;
    g[4] += Bi_4_2 * u_2 + u_0 * Bi_4_0;
    g[5] += u_2 * Bi_5_2 + u_1 * Bi_5_1;
  }  // end of updateGradients

  inline void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_2_2 = dN(ni, 2);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto Bi_4_0 = dN(ni, 2) * icste;
    const auto Bi_4_2 = dN(ni, 0) * icste;
    const auto Bi_5_1 = dN(ni, 2) * icste;
    const auto Bi_5_2 = dN(ni, 1) * icste;
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    Fe[ni_0] += w * (s[3] * Bi_3_0 + s[4] * Bi_4_0 + s[0] * Bi_0_0);
    Fe[ni_1] += w * (s[3] * Bi_3_1 + s[1] * Bi_1_1 + Bi_5_1 * s[5]);
    Fe[ni_2] += w * (Bi_4_2 * s[4] + s[5] * Bi_5_2 + s[2] * Bi_2_2);
  }  // end of updateInnerForces

  inline void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_2_2 = dN(ni, 2);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto Bi_4_0 = dN(ni, 2) * icste;
    const auto Bi_4_2 = dN(ni, 0) * icste;
    const auto Bi_5_1 = dN(ni, 2) * icste;
    const auto Bi_5_2 = dN(ni, 1) * icste;
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto Bj_0_0 = dN(nj, 0);
      const auto Bj_1_1 = dN(nj, 1);
      const auto Bj_2_2 = dN(nj, 2);
      const auto Bj_3_0 = dN(nj, 1) * icste;
      const auto Bj_3_1 = dN(nj, 0) * icste;
      const auto Bj_4_0 = dN(nj, 2) * icste;
      const auto Bj_4_2 = dN(nj, 0) * icste;
      const auto Bj_5_1 = dN(nj, 2) * icste;
      const auto Bj_5_2 = dN(nj, 1) * icste;
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      const auto nj_2 = nj + 2 * nnodes;
      Ke(ni_0, nj_0) +=
          w * (Kip[3] * Bi_0_0 * Bj_3_0 + Kip[27] * Bi_4_0 * Bj_3_0 +
               Bj_0_0 * Kip[18] * Bi_3_0 + Bi_4_0 * Kip[28] * Bj_4_0 +
               Bj_0_0 * Kip[24] * Bi_4_0 + Kip[22] * Bj_4_0 * Bi_3_0 +
               Bj_4_0 * Bi_0_0 * Kip[4] + Kip[21] * Bi_3_0 * Bj_3_0 +
               Bj_0_0 * Kip[0] * Bi_0_0);
      Ke(ni_0, nj_1) +=
          w * (Bj_1_1 * Kip[19] * Bi_3_0 + Bj_3_1 * Kip[3] * Bi_0_0 +
               Bj_1_1 * Bi_4_0 * Kip[25] + Kip[23] * Bj_5_1 * Bi_3_0 +
               Kip[5] * Bj_5_1 * Bi_0_0 + Kip[21] * Bj_3_1 * Bi_3_0 +
               Bj_5_1 * Bi_4_0 * Kip[29] + Kip[27] * Bi_4_0 * Bj_3_1 +
               Kip[1] * Bj_1_1 * Bi_0_0);
      Ke(ni_0, nj_2) +=
          w * (Bj_4_2 * Kip[22] * Bi_3_0 + Bj_4_2 * Bi_0_0 * Kip[4] +
               Bi_4_0 * Bj_5_2 * Kip[29] + Kip[20] * Bj_2_2 * Bi_3_0 +
               Kip[23] * Bj_5_2 * Bi_3_0 + Bj_4_2 * Bi_4_0 * Kip[28] +
               Kip[5] * Bj_5_2 * Bi_0_0 + Kip[2] * Bj_2_2 * Bi_0_0 +
               Bi_4_0 * Bj_2_2 * Kip[26]);
      Ke(ni_1, nj_0) +=
          w * (Kip[10] * Bi_1_1 * Bj_4_0 + Bj_0_0 * Bi_1_1 * Kip[6] +
               Bj_0_0 * Bi_3_1 * Kip[18] + Bi_3_1 * Kip[21] * Bj_3_0 +
               Bi_3_1 * Kip[22] * Bj_4_0 + Bi_1_1 * Kip[9] * Bj_3_0 +
               Bi_5_1 * Kip[33] * Bj_3_0 + Bi_5_1 * Kip[34] * Bj_4_0 +
               Bj_0_0 * Bi_5_1 * Kip[30]);
      Ke(ni_1, nj_1) +=
          w * (Bi_5_1 * Kip[33] * Bj_3_1 + Bi_1_1 * Kip[11] * Bj_5_1 +
               Bi_1_1 * Bj_1_1 * Kip[7] + Kip[23] * Bi_3_1 * Bj_5_1 +
               Bi_3_1 * Kip[21] * Bj_3_1 + Bi_5_1 * Bj_1_1 * Kip[31] +
               Bi_3_1 * Bj_1_1 * Kip[19] + Bi_1_1 * Bj_3_1 * Kip[9] +
               Bi_5_1 * Bj_5_1 * Kip[35]);
      Ke(ni_1, nj_2) +=
          w * (Kip[23] * Bi_3_1 * Bj_5_2 + Bi_5_1 * Kip[35] * Bj_5_2 +
               Bi_1_1 * Kip[11] * Bj_5_2 + Kip[20] * Bi_3_1 * Bj_2_2 +
               Bi_1_1 * Kip[8] * Bj_2_2 + Bi_5_1 * Bj_2_2 * Kip[32] +
               Kip[10] * Bj_4_2 * Bi_1_1 + Bj_4_2 * Bi_3_1 * Kip[22] +
               Bi_5_1 * Bj_4_2 * Kip[34]);
      Ke(ni_2, nj_0) +=
          w * (Bj_0_0 * Kip[24] * Bi_4_2 + Kip[33] * Bi_5_2 * Bj_3_0 +
               Kip[28] * Bi_4_2 * Bj_4_0 + Bi_5_2 * Kip[34] * Bj_4_0 +
               Bj_0_0 * Bi_2_2 * Kip[12] + Bi_2_2 * Bj_4_0 * Kip[16] +
               Bi_2_2 * Kip[15] * Bj_3_0 + Kip[27] * Bi_4_2 * Bj_3_0 +
               Bj_0_0 * Kip[30] * Bi_5_2);
      Ke(ni_2, nj_1) +=
          w * (Bj_5_1 * Bi_4_2 * Kip[29] + Bj_1_1 * Bi_5_2 * Kip[31] +
               Kip[17] * Bj_5_1 * Bi_2_2 + Bj_3_1 * Bi_2_2 * Kip[15] +
               Bj_1_1 * Bi_2_2 * Kip[13] + Bi_5_2 * Bj_5_1 * Kip[35] +
               Kip[27] * Bj_3_1 * Bi_4_2 + Bj_1_1 * Kip[25] * Bi_4_2 +
               Kip[33] * Bi_5_2 * Bj_3_1);
      Ke(ni_2, nj_2) +=
          w * (Kip[14] * Bi_2_2 * Bj_2_2 + Bi_5_2 * Bj_2_2 * Kip[32] +
               Bj_4_2 * Bi_5_2 * Kip[34] + Kip[17] * Bi_2_2 * Bj_5_2 +
               Bj_5_2 * Bi_4_2 * Kip[29] + Bj_4_2 * Bi_2_2 * Kip[16] +
               Bj_4_2 * Kip[28] * Bi_4_2 + Bi_5_2 * Kip[35] * Bj_5_2 +
               Bj_2_2 * Bi_4_2 * Kip[26]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      mgis::raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator

  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateGradients(mgis::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 6>
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 6> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  void
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

  const mfem::IntegrationRule &
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans) {
    const auto order = 2 * Trans.OrderGrad(&el);
    return mfem::IntRules.Get(el.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return getIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator::
      ~OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
