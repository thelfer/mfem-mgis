#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/TransientHeatTransferBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  TransientHeatTransferBehaviourIntegrator::selectIntegrationRule(
      const mfem::FiniteElement &e, const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  TransientHeatTransferBehaviourIntegrator::buildQuadratureSpace(
      const FiniteElementDiscretization &fed, const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  TransientHeatTransferBehaviourIntegrator::
      TransientHeatTransferBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            TransientHeatTransferBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // TransientHeatTransferBehaviourIntegrator

  real TransientHeatTransferBehaviourIntegrator::getIntegrationPointWeight(
      mfem::ElementTransformation &tr,
      const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  TransientHeatTransferBehaviourIntegrator::getIntegrationRule(
      const mfem::FiniteElement &e,
      const mfem::ElementTransformation &t) const {
    return TransientHeatTransferBehaviourIntegrator::selectIntegrationRule(e,
                                                                           t);
  }
  TransientHeatTransferBehaviourIntegrator::RotationMatrix
  TransientHeatTransferBehaviourIntegrator::getRotationMatrix(
      const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  void TransientHeatTransferBehaviourIntegrator::rotateGradients(
      std::span<real>, const RotationMatrix &) {}  // end of rotateGradients

  std::span<const real>
  TransientHeatTransferBehaviourIntegrator::rotateThermodynamicForces(
      std::span<const real> s, const RotationMatrix &) {
    return s;
  }

  void TransientHeatTransferBehaviourIntegrator::rotateTangentOperatorBlocks(
      std::span<real>, const RotationMatrix &) {}

  inline void TransientHeatTransferBehaviourIntegrator::updateGradients(
      std::span<real> &g,
      const mfem::Vector &T,
      const mfem::Vector &N,
      const size_type ni) noexcept {
    g[0] += T[ni] * N[0];
  }  // end of updateGradients

  inline void TransientHeatTransferBehaviourIntegrator::updateInnerForces(
      mfem::Vector &Fe,
      const std::span<const real> &s,
      const mfem::Vector &N,
      const real w,
      const size_type ni) const noexcept {
    Fe[ni] += w * s[0] * N[ni];
  }  // end of updateInnerForces

  inline void TransientHeatTransferBehaviourIntegrator::updateStiffnessMatrix(
      mfem::DenseMatrix &Ke,
      const std::span<const real> &Kip,
      const mfem::Vector &N,
      const real w,
      const size_type ni) const noexcept {
    const auto Ni = N[ni];
    const auto nnodes = N.Size();
    for (size_type nj = 0; nj != nnodes; ++nj) {
      const auto Nj = N[nj];
      Ke(ni, nj) += w * Kip[0] * Ni * Nj;
    }
  }  // end of updateStiffnessMatrix

  bool TransientHeatTransferBehaviourIntegrator::integrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void TransientHeatTransferBehaviourIntegrator::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void TransientHeatTransferBehaviourIntegrator::updateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void TransientHeatTransferBehaviourIntegrator::computeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  TransientHeatTransferBehaviourIntegrator::
      ~TransientHeatTransferBehaviourIntegrator() = default;

}  // end of namespace mfem_mgis
