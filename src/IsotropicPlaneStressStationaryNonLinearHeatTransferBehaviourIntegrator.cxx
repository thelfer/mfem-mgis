#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator

  real IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationPointWeight(mfem::ElementTransformation &tr,
                                const mfem::IntegrationPoint &ip) const
      noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  void
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::setup(
      const real t, const real dt) {
    BehaviourIntegratorBase::setup(t, dt);
    const auto &pev = this->s1.external_state_variables.find("Temperature");
    if (pev == this->s1.external_state_variables.end()) {
      raise(
          "IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrat"
          "or::setup: "
          "external state variable 'Temperature' is not defined");
    }
    if (std::holds_alternative<mgis::span<real>>(pev->second)) {
      this->uesv = std::get<mgis::span<real>>(pev->second).data();
    } else if (std::holds_alternative<std::vector<real>>(pev->second)) {
      this->uesv = std::get<std::vector<real>>(pev->second).data();
    } else {
      raise(
          "IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrat"
          "or::setup: "
          "external state variable 'Temperature' shall not be uniform");
    }
  }  // end of setup

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &u,
                                                     const mfem::Vector &N,
                                                     const size_type o) {
    auto &ev = this->uesv[o];
    ev = 0;
    for (size_type i = 0; i != u.Size(); ++i) {
      ev += u[i] * N[i];
    }
  }

  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      RotationMatrix
      IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
          getRotationMatrix(const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateGradients(mgis::span<real>, const RotationMatrix &) {
  }  // end of rotateGradients

  mgis::span<const real>
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real>, const RotationMatrix &) {}

  inline void
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
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
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    Fe[ni] += w * (dNi_0 * s[0] + s[1] * dNi_1);
  }  // end of updateInnerForces

  inline void
  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
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
      Ke(ni, nj) += w * (Kip[3] * dNj_1 * dNi_1 + Kip[1] * dNi_0 * dNj_1 +
                         Kip[0] * dNj_0 * dNi_0 + dNj_0 * Kip[2] * dNi_1 +
                         (Kip[5] * dNi_1 + Kip[4] * dNi_0) * N[nj]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator::
      ~IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
