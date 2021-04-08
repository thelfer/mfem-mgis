#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator

  real
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationPointWeight(mfem::ElementTransformation &tr,
                                const mfem::IntegrationPoint &ip) const
      noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      setup(const real t, const real dt) {
    BehaviourIntegratorBase::setup(t, dt);
    const auto &pev = this->s1.external_state_variables.find("Temperature");
    if (pev == this->s1.external_state_variables.end()) {
      raise(
          "OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourInt"
          "egrator::setup: "
          "external state variable 'Temperature' is not defined");
    }
    if (mgis::holds_alternative<mgis::span<real>>(pev->second)) {
      this->uesv = mgis::get<mgis::span<real>>(pev->second).data();
    } else if (mgis::holds_alternative<std::vector<real>>(pev->second)) {
      this->uesv = mgis::get<std::vector<real>>(pev->second).data();
    } else {
      raise(
          "OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourInt"
          "egrator::setup: "
          "external state variable 'Temperature' shall not be uniform");
    }
  }  // end of setup

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateExternalStateVariablesFromUnknownsValues(const mfem::Vector &u,
                                                     const mfem::Vector &N,
                                                     const size_type o) {
    auto &ev = this->uesv[o];
    ev = 0;
    for (size_type i = 0; i != u.Size(); ++i) {
      ev += u[i] * N[i];
    }
  }

  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      RotationMatrix
      OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateGradients(mgis::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 3>
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateThermodynamicForces(mgis::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 3> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      rotateTangentOperatorBlocks(mgis::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  inline void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateGradients(mgis::span<real> &g,
                      const mfem::Vector &u,
                      const mfem::DenseMatrix &dN,
                      const size_type ni) noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    const auto u_0 = u[ni];
    g[0] += u_0 * dNi_0;
    g[1] += dNi_1 * u_0;
    g[2] += dNi_2 * u_0;
  }  // end of updateGradients

  inline void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateInnerForces(mfem::Vector &Fe,
                        const mgis::span<const real> &s,
                        const mfem::DenseMatrix &dN,
                        const real w,
                        const size_type ni) const noexcept {
    const auto dNi_0 = dN(ni, 0);
    const auto dNi_1 = dN(ni, 1);
    const auto dNi_2 = dN(ni, 2);
    Fe[ni] += w * (s[1] * dNi_1 + s[2] * dNi_2 + s[0] * dNi_0);
  }  // end of updateInnerForces

  inline void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateStiffnessMatrix(mfem::DenseMatrix &Ke,
                            const mgis::span<const real> &Kip,
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
          w * (dNi_0 * dNj_1 * Kip[1] + dNi_1 * dNj_0 * Kip[3] +
               dNi_1 * Kip[5] * dNj_2 + dNj_1 * Kip[7] * dNi_2 +
               Kip[8] * dNj_2 * dNi_2 + dNi_2 * dNj_0 * Kip[6] +
               dNi_0 * Kip[2] * dNj_2 + dNj_1 * Kip[4] * dNi_1 +
               dNi_0 * dNj_0 * Kip[0] +
               (dNi_0 * Kip[9] + Kip[10] * dNi_1 + Kip[11] * dNi_2) * N[nj]);
    }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
  }    // end of updateStiffnessMatrix

  bool
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void
  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator::
      ~OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
