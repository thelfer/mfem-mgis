#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator

  real
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      setup(const real t, const real dt) {
    BehaviourIntegratorBase::setup(t, dt);
    const auto &pev = this->s1.external_state_variables.find("Temperature");
    if (pev == this->s1.external_state_variables.end()) {
      raise(
          "OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegr"
          "ator::setup: "
          "external state variable 'Temperature' is not defined");
    }
    auto& value = pev->second.value;
    if (std::holds_alternative<std::span<real>>(value)) {
      this->uesv = std::get<std::span<real>>(value).data();
    } else if (std::holds_alternative<std::vector<real>>(value)) {
      this->uesv = std::get<std::vector<real>>(value).data();
    } else {
      raise(
          "OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegr"
          "ator::setup: "
          "external state variable 'Temperature' shall not be uniform");
    }
  }  // end of setup

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &u,
                                                     const mfem::Vector &N,
                                                     const size_type o) {
    auto &ev = this->uesv[o];
    ev = 0;
    for (size_type i = 0; i != u.Size(); ++i) {
      ev += u[i] * N[i];
    }
  }

  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      RotationMatrix
      OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateGradients(std::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 2>
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 2> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  inline void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateGradients(std::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto u_0 = u[ni];
    g[0] += u_0 * dNi_0;
    g[1] += u_0 * dNi_1;
  }  // end of updateGradients

  inline void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const std::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    Fe[ni] += w * (s[1] * dNi_1 + s[0] * dNi_0);
  }  // end of updateInnerForces

  inline void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::Vector &N,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      Ke(ni, nj) += w * (Kip[1] * dNi_0 * dNj_1 + dNi_1 * Kip[3] * dNj_1 +
                         Kip[0] * dNj_0 * dNi_0 + Kip[2] * dNi_1 * dNj_0 +
                         (Kip[4] * dNi_0 + dNi_1 * Kip[5]) * N[nj]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void
  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      ~OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
