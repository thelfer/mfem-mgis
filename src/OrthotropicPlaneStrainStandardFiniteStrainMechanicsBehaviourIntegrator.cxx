#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of
     // OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator

  real OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }
  const mfem::IntegrationRule &
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      RotationMatrix
      OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
          getRotationMatrix(const size_type i) const {
    return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
  }  // end of getRotationMatrix

  void OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateGradients(std::span<real> g, const RotationMatrix &r) {
    this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
  }  // end of rotateGradients

  std::array<real, 5>
  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateThermodynamicForces(std::span<const real> s,
                                const RotationMatrix &r) {
    std::array<real, 5> rs;
    std::copy(s.begin(), s.end(), rs.begin());
    this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
    return rs;
  }

  void OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      rotateTangentOperatorBlocks(std::span<real> Kip,
                                  const RotationMatrix &r) {
    this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                               r.data());
  }

  bool OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      integrate(const mfem::FiniteElement &e,
                mfem::ElementTransformation &tr,
                const mfem::Vector &u,
                const IntegrationType it) {
    return this->implementIntegrate(e, tr, u, it);
  }  // end of integrate

  void OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateResidual(mfem::Vector &Fe,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &u) {
    this->implementUpdateResidual(Fe, e, tr, u);
  }  // end of updateResidual

  void OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      updateJacobian(mfem::DenseMatrix &Ke,
                     const mfem::FiniteElement &e,
                     mfem::ElementTransformation &tr,
                     const mfem::Vector &) {
    this->implementUpdateJacobian(Ke, e, tr);
  }  // end of updateJacobian

  void OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of computeInnerForces

  OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator::
      ~OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis
