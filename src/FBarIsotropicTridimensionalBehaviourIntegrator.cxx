#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/FBarIsotropicTridimensionalBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  FBarIsotropicTridimensionalBehaviourIntegrator::selectIntegrationRule(
      const mfem::FiniteElement &e, const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  FBarIsotropicTridimensionalBehaviourIntegrator::buildQuadratureSpace(
      const FiniteElementDiscretization &fed, const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  FBarIsotropicTridimensionalBehaviourIntegrator::
      FBarIsotropicTridimensionalBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : FBarBehaviourIntegratorCRTPBase<
            FBarIsotropicTridimensionalBehaviourIntegrator,
            Hypothesis::TRIDIMENSIONAL>(buildQuadratureSpace(fed, m),
                                        std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // FBarIsotropicTridimensionalBehaviourIntegrator

  real
  FBarIsotropicTridimensionalBehaviourIntegrator::getIntegrationPointWeight(
      mfem::ElementTransformation &tr,
      const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  FBarIsotropicTridimensionalBehaviourIntegrator::getIntegrationRule(
      const mfem::FiniteElement &e,
      const mfem::ElementTransformation &t) const {
    return FBarIsotropicTridimensionalBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }

  FBarIsotropicTridimensionalBehaviourIntegrator::RotationMatrix
  FBarIsotropicTridimensionalBehaviourIntegrator::getRotationMatrix(
      const size_type) const {
    return RotationMatrix{};
  }  // end of getRotationMatrix

  bool FBarIsotropicTridimensionalBehaviourIntegrator::
      requiresCurrentSolutionForJacobianAssembly() const noexcept {
    return true;
  }  // end of requiresCurrentSolutionForJacobianAssembly

  void FBarIsotropicTridimensionalBehaviourIntegrator::rotateGradients(
      std::span<real>, const RotationMatrix &) {}  // end of rotateGradients

  std::span<const real>
  FBarIsotropicTridimensionalBehaviourIntegrator::rotateThermodynamicForces(
      std::span<const real> s, const RotationMatrix &) {
    return s;
  }

  void
  FBarIsotropicTridimensionalBehaviourIntegrator::rotateTangentOperatorBlocks(
      std::span<real>, const RotationMatrix &) {}

  bool FBarIsotropicTridimensionalBehaviourIntegrator::integrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void FBarIsotropicTridimensionalBehaviourIntegrator::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void FBarIsotropicTridimensionalBehaviourIntegrator::updateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    this->implementUpdateJacobian(Ke, e, tr, u);
  }  // end of updateJacobian

  void FBarIsotropicTridimensionalBehaviourIntegrator::computeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  FBarIsotropicTridimensionalBehaviourIntegrator::
      ~FBarIsotropicTridimensionalBehaviourIntegrator() = default;

}  // end of namespace mfem_mgis
