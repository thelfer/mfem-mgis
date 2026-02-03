#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator

  real
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      setup(const real t, const real dt) {
    BehaviourIntegratorBase::setup(t, dt);
    const auto &pev = this->s1.external_state_variables.find("Temperature");
    if (pev == this->s1.external_state_variables.end()) {
      raise(
          "IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourInteg"
          "rator::setup: "
          "external state variable 'Temperature' is not defined");
    }
    auto &value = pev->second.value;
    if (std::holds_alternative<std::span<real>>(value)) {
      this->uesv = std::get<std::span<real>>(value).data();
    } else if (std::holds_alternative<std::vector<real>>(value)) {
      this->uesv = std::get<std::vector<real>>(value).data();
    } else {
      raise(
          "IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourInteg"
          "rator::setup: "
          "external state variable 'Temperature' shall not be uniform");
    }
  }  // end of setup

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &u,
                                                     const mfem::Vector &N,
                                                     const size_type o) {
    auto &ev = this->uesv[o];
    ev = 0;
    for (size_type i = 0; i != u.Size(); ++i) {
      ev += u[i] * N[i];
    }
  }

  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      RotationMatrix
      IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
          getRotationMatrix(const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateGradients(std::span<real>, const RotationMatrix &) {
  }  // end of rotateGradients

  std::span<const real>
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &) {
    return s;
  }

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real>, const RotationMatrix &) {}

  inline void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateGradients(std::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto u_0 = u[ni];
    g[0] += u_0 * dNi_0;
    g[1] += u_0 * dNi_1;
    g[2] += dNi_2 * u_0;
  }  // end of updateGradients

  inline void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const std::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    Fe[ni] += w * (s[0] * dNi_0 + s[2] * dNi_2 + s[1] * dNi_1);
  }  // end of updateInnerForces

  inline void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const std::span<const real> &Kip,
                            const mfem::Vector &N,
                            const mfem::DenseMatrix &dN,
                            const real w,
                            const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto dNj_0 = dN(nj, 0);
      const auto dNj_1 = dN(nj, 1);
      const auto dNj_2 = dN(nj, 2);
      Ke(ni, nj) +=
          w * (dNi_1 * dNj_0 * Kip[3] + dNi_1 * dNj_1 * Kip[4] +
               Kip[2] * dNj_2 * dNi_0 + dNi_2 * dNj_1 * Kip[7] +
               dNj_0 * Kip[0] * dNi_0 + Kip[1] * dNi_0 * dNj_1 +
               dNi_1 * Kip[5] * dNj_2 + Kip[8] * dNj_2 * dNi_2 +
               dNi_2 * dNj_0 * Kip[6] +
               (dNi_1 * Kip[10] + Kip[11] * dNi_2 + Kip[9] * dNi_0) * N[nj]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void
  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      ~IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
