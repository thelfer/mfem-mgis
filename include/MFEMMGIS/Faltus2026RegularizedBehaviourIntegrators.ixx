/*!
 * \file   MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.ixx
 * \brief
 * \author Thomas Helfer
 * \date   19/03/2026
 */

#ifndef LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX
#define LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX

#include <array>
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  template <Hypothesis H>
  Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::
      Faltus2026RegularizedIsotropicBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr,
          const Parameters &params)
      : Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>(
            fed, m, std::move(b_ptr)),
        alpha(get<real>(may_throw, params, "PenalizationCoefficient")) {
    checkParameters(
        may_throw, params,
        std::map<std::string, std::string>{
            {"PenalizationCoefficient",
             "penalization coefficient of the difference between the "
             "deformation gradient at the centroid of the element and the "
             "deformation gradient at integration points"}});
  }  // end of Faltus2026RegularizedIsotropicBehaviourIntegrator

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::requiresCurrentSolutionForResidualAssembly() const noexcept {
    return true;
  }  // end of requiresCurrentSolutionForResidualAssembly

  template <Hypothesis H>
  bool Faltus2026RegularizedIsotropicBehaviourIntegrator<
      H>::requiresCurrentSolutionForJacobianAssembly() const noexcept {
    return false;
  }  // end of requiresCurrentSolutionForJacobianAssembly

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::updateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::updateResidual(
        Fe, e, tr, u);
    constexpr auto gsize = []() constexpr {
      if constexpr ((H == Hypothesis::PLANESTRAIN) ||
                    (H == Hypothesis::PLANESTRESS)) {
        return size_type{5};
      } else {
        return size_type{9};
      }
    }
    ();
#ifdef MFEM_THREAD_SAFE
    mfem::DenseMatrix dshape, dshape0;
    dshape.setSize(e.GetDof(), e.GetDim());
    dshape0.setSize(e.GetDof(), e.GetDim());
#else
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    // this->dshape is updated by
    // Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>
#endif
    const auto nnodes = e.GetDof();
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    auto F0 = std::array<real, gsize>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    {  // computation of the mean value of the deformation gradient
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape0);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        auto F0v = std::span<real>(F0);
        this->updateGradients(F0v, u, dshape0, ni);
      }
    }
    // contribution to the residual
    auto dF = std::array<real, gsize>{};
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the gradients of the shape functions
#ifdef MFEM_THREAD_SAFE
      e.CalcPhysDShape(tr, dshape);
#else
      e.CalcPhysDShape(tr, this->dshape);
#endif
      // offset of the integration point
      const auto o = eoffset + i;
      // current estimation of the deformation gradient at the end of the time
      // step
      auto F = this->s1.gradients.subspan(o * gsize, gsize);
      for (size_type c = 0; c != gsize; ++c) {
        dF[c] = this->alpha * (F[c] - F0[c]);
      }
      const auto w = this->getIntegrationPointWeight(tr, ip);
      for (size_type ni = 0; ni != nnodes; ++ni) {
#ifdef MFEM_THREAD_SAFE
        this->updateInnerForces(Fe, dF, dshape, w, ni);
#else
        this->updateInnerForces(Fe, dF, this->dshape, w, ni);
#endif
        this->updateInnerForces(Fe, dF, dshape0, -w, ni);
      }
    }
  }  // end of updateResidual

  template <Hypothesis H>
  void Faltus2026RegularizedIsotropicBehaviourIntegrator<H>::updateJacobian(
      mfem::DenseMatrix &Je,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>::updateJacobian(
        Je, e, tr, u);
    constexpr auto gsize = []() constexpr {
      if constexpr ((H == Hypothesis::PLANESTRAIN) ||
                    (H == Hypothesis::PLANESTRESS)) {
        return size_type{5};
      } else {
        return size_type{9};
      }
    }
    ();
#ifdef MFEM_THREAD_SAFE
    mfem::DenseMatrix dshape, dshape0;
    dshape.setSize(e.GetDof(), e.GetDim());
    dshape0.setSize(e.GetDof(), e.GetDim());
#else
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    // this->dshape is updated by
    // Faltus2026RegularizedIsotropicBehaviourIntegratorBase<H>
#endif
    {  // computation of the mean value of the deformation gradient
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape0);
    }
    const auto &ir = this->getIntegrationRule(e, tr);
    auto Km = std::array<real, gsize * gsize>{};
    for (size_type i = 0; i != gsize; ++i) {
      for (size_type j = 0; j != gsize; ++j) {
        Km[i * gsize + j] = 0;
      }
      Km[i * gsize + i] = this->alpha;
    }
    const auto nnodes = e.GetDof();
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // get the gradients of the shape functions
#ifdef MFEM_THREAD_SAFE
      e.CalcPhysDShape(tr, dshape);
#else
      e.CalcPhysDShape(tr, this->dshape);
#endif
      const auto w = this->getIntegrationPointWeight(tr, ip);
      for (size_type ni = 0; ni != nnodes; ++ni) {
#ifdef MFEM_THREAD_SAFE
        this->updateStiffnessMatrix(Je, Km, dshape, dshape, w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape0, dshape0, -w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape, dshape0, -w, ni);
        this->updateStiffnessMatrix(Je, Km, dshape0, dshape, w, ni);
#else
        this->updateStiffnessMatrix(Je, Km, this->dshape, this->dshape, w, ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape0, this->dshape0, -w,
                                    ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape, this->dshape0, -w,
                                    ni);
        this->updateStiffnessMatrix(Je, Km, this->dshape0, this->dshape, w, ni);
#endif
      }
    }
  }  // end of updateJacobian

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_FALTUS2026REGULARIZEDBEHAVIOURINTEGRATORS_IXX */
