/*!
 * \file   TridimensionalMicromorphicDamageBehaviourIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   31/03/2023
 */

#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/TridimensionalMicromorphicDamageBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  TridimensionalMicromorphicDamageBehaviourIntegrator::selectIntegrationRule(
      const mfem::FiniteElement &e, const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  TridimensionalMicromorphicDamageBehaviourIntegrator::buildQuadratureSpace(
      const FiniteElementDiscretization &fed, const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  TridimensionalMicromorphicDamageBehaviourIntegrator::
      TridimensionalMicromorphicDamageBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(buildQuadratureSpace(fed, m),
                                std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ISOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of TridimensionalMicromorphicDamageBehaviourIntegrator

  real TridimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationPointWeight(
          mfem::ElementTransformation &tr,
          const mfem::IntegrationPoint &ip) const noexcept {
    return ip.weight * tr.Weight();
  }  // end of getIntegrationPointWeight

  const mfem::IntegrationRule &
  TridimensionalMicromorphicDamageBehaviourIntegrator::getIntegrationRule(
      const mfem::FiniteElement &e,
      const mfem::ElementTransformation &t) const {
    return TridimensionalMicromorphicDamageBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }  // end of getIntegrationRule

  bool TridimensionalMicromorphicDamageBehaviourIntegrator::integrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &d_chi,
      const IntegrationType it) {
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    // integration rule
    const auto &ir = this->getIntegrationRule(e, tr);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    shape.SetSize(e.GetDof());
    mfem::DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    this->shape.SetSize(e.GetDof());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    const auto gsize = this->s1.gradients_stride;
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      const auto o = eoffset + i;
      // get the values of the shape functions
      e.CalcPhysShape(tr, shape);
      // get the gradients of the shape functions
      e.CalcPhysDShape(tr, dshape);
      //
      auto g = this->s1.gradients.subspan(o * gsize, gsize);
      std::copy(this->macroscopic_gradients.begin(),
                this->macroscopic_gradients.end(), g.begin());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // computations of the gradients:
        // - g[0]: micromorphic damage
        // - g[1]: first component of the gradient of the micromorphic damage
        // - g[2]: second component of the gradient of the micromorphic damage
        g[0] += d_chi[ni] * shape[ni];
        g[1] += d_chi[ni] * dshape(ni, 0);
        g[2] += d_chi[ni] * dshape(ni, 1);
        g[3] += d_chi[ni] * dshape(ni, 2);
      }
      if (!this->performsLocalBehaviourIntegration(o, it)) {
        return false;
      }
    }
    return true;
  }  // end of integrate

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    shape.SetSize(e.GetDof());
    mfem::DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    this->shape.SetSize(e.GetDof());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    const auto thsize = this->s1.thermodynamic_forces_stride;
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe.SetSize(e.GetDof());
    Fe = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the values of the shape functions
      e.CalcPhysShape(tr, shape);
      // get the gradients of the shape functions
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = this->getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // s contains:
        // - a_chi, the dual force conjugated with d_chi (scalar)
        // - b_chi, the dual force conjugated with the gradient of d_chi
        // (vector)
        Fe[ni] += w * (s[0] * shape[ni] +      //
                       s[1] * dshape(ni, 0) +  //
                       s[2] * dshape(ni, 1) +  //
                       s[3] * dshape(ni, 2));
      }
    }
  }  // end of updateResidual

  void TridimensionalMicromorphicDamageBehaviourIntegrator::updateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    shape.SetSize(e.GetDof());
    mfem::DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    this->shape.SetSize(e.GetDof());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto &ir = this->getIntegrationRule(e, tr);
    const auto nnodes = e.GetDof();
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Ke.SetSize(e.GetDof(), e.GetDof());
    Ke = 0.;
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // get the gradients of the shape functions
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      e.CalcPhysShape(tr, shape);
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = this->getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      const auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // Kip contains:
        // - the derivative of a_chi with respect to d_chi (scalar)
        // - the derivative of b_chi with respect to the gradient of d_chi (3 *
        // 3 matrix)
        const auto Ni = shape[ni];
        const auto dNi_0 = dshape(ni, 0);
        const auto dNi_1 = dshape(ni, 1);
        const auto dNi_2 = dshape(ni, 2);
        for (size_type nj = 0; nj != nnodes; ++nj) {
          const auto Nj = shape[nj];
          const auto dNj_0 = dshape(nj, 0);
          const auto dNj_1 = dshape(nj, 1);
          const auto dNj_2 = dshape(nj, 2);
          Ke(ni, nj) += w * (Kip[0] * Ni * Nj +  //
                             dNi_1 * dNj_0 * Kip[4] + dNi_1 * dNj_1 * Kip[5] +
                             Kip[3] * dNj_2 * dNi_0 + dNi_2 * dNj_1 * Kip[8] +
                             dNj_0 * Kip[1] * dNi_0 + Kip[2] * dNi_0 * dNj_1 +
                             dNi_1 * Kip[6] * dNj_2 + Kip[9] * dNj_2 * dNi_2 +
                             dNi_2 * dNj_0 * Kip[7]);
        }
      }
    }
  }  // end of updateJacobian

  void TridimensionalMicromorphicDamageBehaviourIntegrator::computeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    shape.SetSize(e.GetDof());
    mfem::DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    this->shape.SetSize(e.GetDof());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    const auto thsize = this->s1.thermodynamic_forces_stride;
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe.SetSize(e.GetDof());
    Fe = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the values of the shape functions
      e.CalcPhysShape(tr, shape);
      // get the gradients of the shape functions
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = this->getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // s contains:
        // - a_chi, the dual force conjugated with d_chi (scalar)
        // - b_chi, the dual force conjugated with the gradient of d_chi
        // (vector)
        Fe[ni] += w * (s[1] * dshape(ni, 0) +  //
                       s[2] * dshape(ni, 1) +  //
                       s[3] * dshape(ni, 2));
      }
    }
  }  // end of computeInnerForces

  TridimensionalMicromorphicDamageBehaviourIntegrator::
      ~TridimensionalMicromorphicDamageBehaviourIntegrator() = default;

}  // end of namespace mfem_mgis