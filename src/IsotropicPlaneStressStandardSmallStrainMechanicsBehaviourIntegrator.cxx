#include "MFEMMGIS/IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  inline void
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    g[0] += Bi_0_0 * u_0;
    g[1] += u_1 * Bi_1_1;
    g[2] += 0;
    g[3] += Bi_3_0 * u_0 + u_1 * Bi_3_1;
  }  // end of updateGradients

  inline void
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
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
    Fe[ni_0] += Bi_0_0 * s[0] * w + s[3] * w * Bi_3_0;
    Fe[ni_1] += Bi_3_1 * s[3] * w + Bi_1_1 * w * s[1];
  }  // end of updateInnerForces

  inline void
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
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
          Kip[0] * Bj_0_0 * w * Bi_0_0 + Bj_0_0 * Bi_3_0 * w * Kip[12] +
          Bi_3_0 * w * Bj_3_0 * Kip[15] + Kip[3] * w * Bi_0_0 * Bj_3_0;
      Ke(ni_0, nj_1) +=
          Bj_3_1 * Kip[3] * w * Bi_0_0 + Bj_3_1 * Bi_3_0 * w * Kip[15] +
          Kip[13] * Bi_3_0 * Bj_1_1 * w + Kip[1] * Bj_1_1 * w * Bi_0_0;
      Ke(ni_1, nj_0) +=
          w * Bi_3_1 * Bj_3_0 * Kip[15] + Bi_1_1 * Bj_0_0 * Kip[4] * w +
          Bi_1_1 * Kip[7] * w * Bj_3_0 + Bj_0_0 * w * Bi_3_1 * Kip[12];
      Ke(ni_1, nj_1) +=
          Bj_3_1 * w * Bi_3_1 * Kip[15] + Bi_1_1 * Bj_3_1 * Kip[7] * w +
          Kip[13] * Bj_1_1 * w * Bi_3_1 + Bi_1_1 * Bj_1_1 * w * Kip[5];
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      mgis::raise("invalid behaviour symmetry");
    }
  }  // end of
     // IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      setRotationMatrix(const RotationMatrix2D &) {
    mgis::raise(
        "IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::"
        "setRotationMatrix: invalid call");
  }

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      setRotationMatrix(const RotationMatrix3D &) {
    mgis::raise(
        "IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::"
        "setRotationMatrix: invalid call");
  }

  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
          getRotationMatrix() const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateGradients(mgis::span<real>, const RotationMatrix &) {
  }  // end of rotateGradients

  mgis::span<const real>
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real>, const RotationMatrix &) {}

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr,
                         const mfem::Vector &u) {
    this->implementComputeInnerForces(Fe, e, tr, u);
  }  // end of computeInnerForces

  void IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      computeStiffnessMatrix(mfem::DenseMatrix &Ke,
                             const mfem::FiniteElement &e,
                             mfem::ElementTransformation &tr,
                             const mfem::Vector &) {
    this->implementComputeStiffnessMatrix(Ke, e, tr);
  }  // end of computeStiffnessMatrix

  const mfem::IntegrationRule &
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &el,
                         const mfem::ElementTransformation &Trans) {
    const auto order = 2 * Trans.OrderGrad(&el);
    return mfem::IntRules.Get(el.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return getIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator::
      ~IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
