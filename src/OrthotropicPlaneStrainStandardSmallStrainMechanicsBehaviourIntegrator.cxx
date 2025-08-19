#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator

  real OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateGradients(std::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 4>
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 4> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  inline void
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      updateGradients(std::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto nnodes = dN.NumRows();
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    g[0] += u_0 * Bi_0_0;
    g[1] += Bi_1_1 * u_1;
    g[2] += 0;
    g[3] += u_1 * Bi_3_1 + u_0 * Bi_3_0;
  }  // end of updateGradients

  inline void
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const std::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto nnodes = dN.NumRows();
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    Fe[ni_0] += w * (Bi_3_0 * s[3] + s[0] * Bi_0_0);
    Fe[ni_1] += w * (Bi_3_1 * s[3] + Bi_1_1 * s[1]);
  }  // end of updateInnerForces

  inline void
  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto Bj_0_0 = dN(nj, 0);
      const auto Bj_1_1 = dN(nj, 1);
      const auto Bj_3_0 = dN(nj, 1) * icste;
      const auto Bj_3_1 = dN(nj, 0) * icste;
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      Ke(ni_0, nj_0) +=
          w * (Kip[12] * Bj_0_0 * Bi_3_0 + Bi_0_0 * Kip[0] * Bj_0_0 +
               Bi_0_0 * Bj_3_0 * Kip[3] + Bj_3_0 * Kip[15] * Bi_3_0);
      Ke(ni_0, nj_1) +=
          w * (Kip[15] * Bj_3_1 * Bi_3_0 + Kip[13] * Bi_3_0 * Bj_1_1 +
               Bi_0_0 * Kip[3] * Bj_3_1 + Bi_0_0 * Kip[1] * Bj_1_1);
      Ke(ni_1, nj_0) +=
          w * (Bi_1_1 * Bj_0_0 * Kip[4] + Bi_3_1 * Kip[12] * Bj_0_0 +
               Bj_3_0 * Bi_1_1 * Kip[7] + Bi_3_1 * Bj_3_0 * Kip[15]);
      Ke(ni_1, nj_1) +=
          w * (Bi_3_1 * Kip[15] * Bj_3_1 + Bi_3_1 * Kip[13] * Bj_1_1 +
               Bi_1_1 * Bj_3_1 * Kip[7] + Kip[5] * Bi_1_1 * Bj_1_1);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator::
      ~OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
