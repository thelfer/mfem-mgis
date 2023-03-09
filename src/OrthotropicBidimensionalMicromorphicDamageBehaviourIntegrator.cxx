/*!
 * \file   OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   07/12/2021
 */

#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.hxx"

namespace mfem_mgis {

  const mfem::IntegrationRule &
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      selectIntegrationRule(const mfem::FiniteElement &e,
                            const mfem::ElementTransformation &t) {
    const auto order = 2 * t.OrderGrad(&e);
    return mfem::IntRules.Get(e.GetGeomType(), order);
  }

  std::shared_ptr<const PartialQuadratureSpace>
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      buildQuadratureSpace(const FiniteElementDiscretization &fed,
                           const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return selectIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(buildQuadratureSpace(fed, m),
                                std::move(b_ptr)) {
    if (this->b.symmetry != Behaviour::ORTHOTROPIC) {
      raise("invalid behaviour symmetry");
    }
  }  // end of OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator

  real OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationPointWeight(mfem::ElementTransformation &tr,
                                const mfem::IntegrationPoint &ip) const
      noexcept {
    return ip.weight * tr.Weight();
  }  // end of getIntegrationPointWeight

  const mfem::IntegrationRule &
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      getIntegrationRule(const mfem::FiniteElement &e,
                         const mfem::ElementTransformation &t) const {
    return OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
        selectIntegrationRule(e, t);
  }  // end of getIntegrationRule

  bool OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::integrate(
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
      }
      const auto r = this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
      this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());
      if (!this->performsLocalBehaviourIntegration(o, it)) {
        return false;
      }
    }
    return true;
  }  // end of integrate

  void
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::updateResidual(
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
      const auto r = this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
      std::array<real, 2> rs;
      std::copy(s.begin(), s.end(), rs.begin());
      this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // s contains:
        // - a_chi, the dual force conjugated with d_chi (scalar)
        // - b_chi, the dual force conjugated with the gradient of d_chi
        // (vector)
        Fe[ni] += w * (rs[0] * shape[ni] +      //
                       rs[1] * dshape(ni, 0) +  //
                       rs[2] * dshape(ni, 1));
      }
    }
  }  // end of updateResidual

  void
  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::updateJacobian(
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
      const auto r = this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
      this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(),
                                                 r.data());
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // Kip contains:
        // - the derivative of a_chi with respect to d_chi (scalar)
        // - the derivative of b_chi with respect to the gradient of d_chi (2 *
        // 2 matrix)
        const auto Ni = shape[ni];
        const auto dNi_0 = dshape(ni, 0);
        const auto dNi_1 = dshape(ni, 1);
        for (size_type nj = 0; nj != nnodes; ++nj) {
          const auto Nj = shape[nj];
          const auto dNj_0 = dshape(nj, 0);
          const auto dNj_1 = dshape(nj, 1);
          Ke(ni, nj) += w * (Kip[0] * Ni * Nj +        //
                             Kip[2] * dNi_0 * dNj_1 +  //
                             dNj_0 * Kip[3] * dNi_1 +  //
                             dNj_1 * dNi_1 * Kip[4] +  //
                             dNj_0 * dNi_0 * Kip[1]);
        }
      }
    }
  }  // end of updateJacobian

  void OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      computeInnerForces(mfem::Vector &Fe,
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
      const auto r = this->get_rotation_fct_ptr(this->r2D, this->r3D, i);
      std::array<real, 2> rs;
      std::copy(s.begin(), s.end(), rs.begin());
      this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), r.data());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        // s contains:
        // - a_chi, the dual force conjugated with d_chi (scalar)
        // - b_chi, the dual force conjugated with the gradient of d_chi
        // (vector)
        Fe[ni] += w * (rs[1] * dshape(ni, 0) +  //
                       rs[2] * dshape(ni, 1));
      }
    }
  }  // end of computeInnerForces

  OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator::
      ~OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator() =
          default;

}  // end of namespace mfem_mgis