#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator

  real
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationPointWeight(mfem::ElementTransformation &tr,
                                const mfem::IntegrationPoint &ip) const
      noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateGradients(mgis::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 9>
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 9> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  inline void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto nnodes = dN.NumRows();
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    const auto u_2 = u[ni + 2 * nnodes];
    g[0] += u_0 * dNi_0;
    g[1] += u_1 * dNi_1;
    g[2] += u_2 * dNi_2;
    g[3] += u_0 * dNi_1;
    g[4] += dNi_0 * u_1;
    g[5] += u_0 * dNi_2;
    g[6] += dNi_0 * u_2;
    g[7] += u_1 * dNi_2;
    g[8] += dNi_1 * u_2;
  }  // end of updateGradients

  inline void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto nnodes = dN.NumRows();
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    Fe[ni_0] += w * (s[5] * dNi_2 + dNi_1 * s[3] + s[0] * dNi_0);
    Fe[ni_1] += w * (dNi_0 * s[4] + s[1] * dNi_1 + s[7] * dNi_2);
    Fe[ni_2] += w * (dNi_1 * s[8] + s[2] * dNi_2 + dNi_0 * s[6]);
  }  // end of updateInnerForces

  inline void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      const auto dNj_2 = dN(nj, 2);
      const auto nj_0 = nj;
      const auto nj_1 = nj + nnodes;
      const auto nj_2 = nj + 2 * nnodes;
      Ke(ni_0, nj_0) += w * (dNj_2 * dNi_2 * Kip[50] + dNj_0 * Kip[0] * dNi_0 +
                             dNj_2 * Kip[5] * dNi_0 + Kip[48] * dNi_2 * dNj_1 +
                             Kip[45] * dNi_2 * dNj_0 + dNj_1 * dNi_0 * Kip[3] +
                             dNi_1 * Kip[30] * dNj_1 + dNi_1 * Kip[27] * dNj_0 +
                             dNj_2 * dNi_1 * Kip[32]);
      Ke(ni_0, nj_1) += w * (dNj_2 * dNi_1 * Kip[34] + dNi_1 * Kip[28] * dNj_1 +
                             Kip[1] * dNj_1 * dNi_0 + dNj_2 * Kip[52] * dNi_2 +
                             Kip[49] * dNi_2 * dNj_0 + Kip[46] * dNi_2 * dNj_1 +
                             dNj_2 * Kip[7] * dNi_0 + Kip[4] * dNj_0 * dNi_0 +
                             dNi_1 * dNj_0 * Kip[31]);
      Ke(ni_0, nj_2) += w * (Kip[35] * dNi_1 * dNj_1 + dNi_1 * Kip[33] * dNj_0 +
                             Kip[8] * dNj_1 * dNi_0 + dNj_2 * dNi_2 * Kip[47] +
                             dNi_2 * Kip[53] * dNj_1 + dNj_0 * dNi_0 * Kip[6] +
                             dNj_2 * Kip[2] * dNi_0 + dNi_2 * dNj_0 * Kip[51] +
                             dNj_2 * dNi_1 * Kip[29]);
      Ke(ni_1, nj_0) += w * (dNj_2 * dNi_0 * Kip[41] + dNi_2 * dNj_0 * Kip[63] +
                             dNi_2 * Kip[66] * dNj_1 + dNi_1 * dNj_0 * Kip[9] +
                             dNi_1 * Kip[12] * dNj_1 + dNj_2 * Kip[68] * dNi_2 +
                             dNj_2 * dNi_1 * Kip[14] + Kip[39] * dNj_1 * dNi_0 +
                             Kip[36] * dNj_0 * dNi_0);
      Ke(ni_1, nj_1) += w * (dNj_2 * Kip[43] * dNi_0 + Kip[13] * dNi_1 * dNj_0 +
                             Kip[64] * dNi_2 * dNj_1 + Kip[67] * dNi_2 * dNj_0 +
                             Kip[70] * dNj_2 * dNi_2 + dNj_0 * Kip[40] * dNi_0 +
                             Kip[37] * dNj_1 * dNi_0 + Kip[10] * dNi_1 * dNj_1 +
                             dNj_2 * dNi_1 * Kip[16]);
      Ke(ni_1, nj_2) += w * (Kip[71] * dNi_2 * dNj_1 + dNj_2 * dNi_0 * Kip[38] +
                             dNj_2 * Kip[65] * dNi_2 + dNi_1 * dNj_0 * Kip[15] +
                             Kip[42] * dNj_0 * dNi_0 + dNi_1 * Kip[17] * dNj_1 +
                             dNj_2 * dNi_1 * Kip[11] + dNi_2 * dNj_0 * Kip[69] +
                             dNj_1 * Kip[44] * dNi_0);
      Ke(ni_2, nj_0) += w * (dNj_0 * dNi_0 * Kip[54] + dNj_2 * Kip[23] * dNi_2 +
                             dNi_2 * dNj_0 * Kip[18] + dNj_2 * dNi_1 * Kip[77] +
                             dNi_1 * dNj_0 * Kip[72] + dNi_2 * Kip[21] * dNj_1 +
                             dNj_1 * Kip[57] * dNi_0 + dNi_1 * Kip[75] * dNj_1 +
                             dNj_2 * Kip[59] * dNi_0);
      Ke(ni_2, nj_1) += w * (dNi_1 * dNj_0 * Kip[76] + Kip[58] * dNj_0 * dNi_0 +
                             Kip[61] * dNj_2 * dNi_0 + dNi_2 * dNj_1 * Kip[19] +
                             Kip[55] * dNj_1 * dNi_0 + dNj_2 * dNi_2 * Kip[25] +
                             Kip[73] * dNi_1 * dNj_1 + dNj_2 * dNi_1 * Kip[79] +
                             dNi_2 * dNj_0 * Kip[22]);
      Ke(ni_2, nj_2) += w * (dNj_2 * Kip[56] * dNi_0 + dNj_0 * Kip[60] * dNi_0 +
                             dNi_1 * Kip[80] * dNj_1 + Kip[24] * dNi_2 * dNj_0 +
                             Kip[26] * dNi_2 * dNj_1 + dNj_2 * dNi_1 * Kip[74] +
                             dNi_1 * dNj_0 * Kip[78] + dNj_2 * Kip[20] * dNi_2 +
                             Kip[62] * dNj_1 * dNi_0);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void
  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator::
      ~OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
