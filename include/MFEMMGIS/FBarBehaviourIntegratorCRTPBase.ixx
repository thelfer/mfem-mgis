/*!
 * \file   include/MFEMMGIS/FBarBehaviourIntegratorCRTPBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX
#define LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX

#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Raise.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"

namespace mfem_mgis {

  template <typename Child, Hypothesis H>
  FBarBehaviourIntegratorCRTPBase<Child, H>::FBarBehaviourIntegratorCRTPBase(
      std::shared_ptr<const PartialQuadratureSpace> s,
      std::unique_ptr<const Behaviour> b_ptr)
      : BehaviourIntegratorBase(s, std::move(b_ptr)) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr auto usize = static_cast<size_type>(Traits::unknownsSize);
    const auto &fed = s->getFiniteElementDiscretization();
    const auto nc = getNumberOfComponents(fed);
    if (nc != usize) {
      raise("the number of components of the finite element space (" +
            std::to_string(nc) +
            ") does not match "
            "the number of components expected by the behaviour integrator (" +
            std::to_string(usize) + ")");
    }
  }  // end of FBarBehaviourIntegratorCRTPBase

  template <typename Child, Hypothesis H>
  bool FBarBehaviourIntegratorCRTPBase<Child, H>::implementIntegrate(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u,
      const IntegrationType it) {
    //
    using Traits = BehaviourIntegratorTraits<Child>;
    static_assert(
        ((Traits::gradientsComputationRequiresShapeFunctions) ||
         (Traits::gradientsComputationRequiresShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shape functions are required to "
                  "compute the gradients");
    constexpr auto dimension = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    // inverse of the exponent of the FBar regularisation
    constexpr auto ie = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    constexpr auto gsize = tfel::math::TensorDimeToSize<dimension>::value;
    //
    auto &child = static_cast<Child &>(*this);
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    const auto &ir = child.getIntegrationRule(e, tr);
    if ((it == IntegrationType::PREDICTION_TANGENT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_SECANT_OPERATOR) ||
        (it == IntegrationType::PREDICTION_ELASTIC_OPERATOR)) {
      for (size_type i = 0; i != ir.GetNPoints(); ++i) {
        // offset of the integration point
        const auto o = eoffset + i;
        if (!this->performsLocalBehaviourIntegration(o, it)) {
          return false;
        }
        // rotate the tangent operator blocks
        const auto r = child.getRotationMatrix(o);
        auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
        child.rotateTangentOperatorBlocks(Kip, r);
      }
      return true;
    }
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape0;
    mfem::DenseMatrix dshape0;
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape0.SetSize(e.GetDof());
      shape.SetSize(e.GetDof());
    }
    dshape0.setSize(e.GetDof(), e.GetDim());
    dshape.setSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape0.SetSize(e.GetDof());
      this->shape.SetSize(e.GetDof());
    }
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    auto F0 = tfel::math::tensor<dimension, real>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    {  // computation of the mean value of the deformation gradient
      auto F0v = std::span<real>(F0.data(), F0.size());
      const auto &ir0 = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir0.IntPoint(0);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        e.CalcPhysShape(tr, shape0);
      }
      e.CalcPhysDShape(tr, dshape0);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        child.updateGradients(F0v, u, dshape0, ni);
      }
    }
    const auto J0 = tfel::math::det(F0);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      // offset of the integration point
      const auto o = eoffset + i;
      if constexpr (evaluateShapeFunctions) {
        // get the gradients of the shape functions
        e.CalcPhysShape(tr, shape);
      }
      if constexpr (Traits::updateExternalStateVariablesFromUnknownsValues) {
        child.updateExternalStateVariablesFromUnknownsValues(u, shape, o);
      }
      // get the gradients of the shape functions
      e.CalcPhysDShape(tr, dshape);
      auto F = tfel::math::tensor<dimension, real>{};
      std::copy(this->macroscopic_gradients.begin(),
                this->macroscopic_gradients.end(), F.begin());
      auto Fv = std::span<real>(F.data(), F.size());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
          child.updateGradients(Fv, u, shape, dshape, ni);
        } else {
          child.updateGradients(Fv, u, dshape, ni);
        }
      }
      const auto J = tfel::math::det(F);
      const auto alpha = tfel::math::power<1, ie>(J0 / J);
      auto g = this->s1.gradients.subspan(o * gsize, gsize);
      if constexpr (H == Hypothesis::PLANESTRAIN) {
        for (size_type c = 0; c != gsize; ++c) {
          g[c] = alpha * F[c];
        }
        g[2] = 1;
      } else {
        for (size_type c = 0; c != gsize; ++c) {
          g[c] = alpha * F[c];
        }
      }
      const auto r = child.getRotationMatrix(o);
      child.rotateGradients(g, r);
      if (!this->performsLocalBehaviourIntegration(o, it)) {
        return false;
      }
      // Here we rotate the tangent operator blocks but not the thermodynamic
      // forces.
      if (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
        auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
        child.rotateTangentOperatorBlocks(Kip, r);
      }
    }
    return true;
  }  // end of implementIntegrate

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementUpdateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of implementUpdateResidual

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementComputeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shap functions are required to "
                  "compute the inner forces");
    auto &child = static_cast<Child &>(*this);
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    dshape.setSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    const auto thsize = this->s1.thermodynamic_forces_stride;
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe.SetSize(e.GetDof() * Traits::unknownsSize);
    Fe = 0.;
    const auto &ir = child.getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        // get the shape functions
        e.CalcPhysShape(tr, shape);
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        // get the gradients of the shape functions
        e.CalcPhysDShape(tr, dshape);
      }
      // get the weights associated to point ip
      const auto w = child.getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      const auto r = child.getRotationMatrix(o);
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      const auto &rs = child.rotateThermodynamicForces(s, r);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (evaluateShapeFunctionsDerivatives) {
          child.updateInnerForces(Fe, rs, dshape, w, ni);
        } else {
          child.updateInnerForces(Fe, rs, shape, w, ni);
        }
      }
    }
  }  // end of implementComputeInnerForces

  template <typename Child, Hypothesis H>
  void FBarBehaviourIntegratorCRTPBase<Child, H>::implementUpdateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &u) {
    using Traits = BehaviourIntegratorTraits<Child>;
    auto &child = static_cast<Child &>(*this);
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    constexpr auto dimension = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    // inverse of the exponent of the FBar regularisation
    constexpr auto ie = []() constexpr->unsigned short {
      return (H == Hypothesis::PLANESTRAIN) ? 2 : 3;
    }
    ();
    constexpr auto gsize = tfel::math::TensorDimeToSize<dimension>::value;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the jacobian, this is not supported");
    static_assert(evaluateShapeFunctionsDerivatives,
                  "the derivatives of the shape functions are required to "
                  "compute the jacobian, this is not supported");
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape0, shape;
    mfem::DenseMatrix dshape0, dshape;
    if constexpr (evaluateShapeFunctions) {
      shape0.SetSize(e.GetDof());
      shape.SetSize(e.GetDof());
    }
    dshape0.SetSize(e.GetDof(), e.GetDim());
    dshape.SetSize(e.GetDof(), e.GetDim());
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape0.SetSize(e.GetDof());
      this->shape.SetSize(e.GetDof());
    }
    this->dshape0.SetSize(e.GetDof(), e.GetDim());
    this->dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    const auto nnodes = e.GetDof();
    auto F0 = tfel::math::tensor<dimension, real>{};
    std::copy(this->macroscopic_gradients.begin(),
              this->macroscopic_gradients.end(), F0.begin());
    {  // computation of the mean value of the deformation gradient
      auto F0v = std::span<real>(F0.data(), F0.size());
      const auto &ir = mfem::IntRules.Get(e.GetGeomType(), 0);
      const auto &ip = ir.IntPoint(0);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        e.CalcPhysShape(tr, shape0);
      }
      e.CalcPhysDShape(tr, dshape0);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        child.updateGradients(F0v, u, dshape0, ni);
      }
    }
    const auto J0 = tfel::math::det(F0);
    // element offset
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Ke.SetSize(e.GetDof() * Traits::unknownsSize,
               e.GetDof() * Traits::unknownsSize);
    Ke = 0.;
    const auto &ir = child.getIntegrationRule(e, tr);
    for (size_type i = 0; i != ir.GetNPoints(); ++i) {
      // get the gradients of the shape functions
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      if constexpr (evaluateShapeFunctions) {
        // get the shape functions
        e.CalcPhysShape(tr, shape);
      }
      // get the derivatives of the shape functions
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = child.getIntegrationPointWeight(tr, ip);
      // offset of the integration point
      const auto o = eoffset + i;
      auto F = tfel::math::tensor<dimension, real>{};
      std::copy(this->macroscopic_gradients.begin(),
                this->macroscopic_gradients.end(), F.begin());
      auto Fv = std::span<real>(F.data(), F.size());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (Traits::gradientsComputationRequiresShapeFunctions) {
          child.updateGradients(Fv, u, shape, dshape, ni);
        } else {
          child.updateGradients(Fv, u, dshape, ni);
        }
      }
      const auto J = tfel::math::det(F);
      const auto alpha = tfel::math::power<1, ie>(J0 / J);
      auto dJ_dF = tfel::math::computeDeterminantDerivative(F);
      auto dJ_dF0 = tfel::math::computeDeterminantDerivative(F0);
      const auto dalpha_dF0 = (alpha / (ie * J0)) * dJ_dF0;
      const auto dalpha_dF = (-alpha / (ie * J)) * dJ_dF;
      auto dFbar_dF0 = tfel::math::eval(F ^ dalpha_dF0);
      auto dFbar_dF = tfel::math::eval(
          alpha * tfel::math::t2tot2<dimension, real>::Id() + (F ^ dalpha_dF));
      if constexpr (ie == 2u) {
        // those are brutal corrections that seems working
        // we shall verify this someday
        for (size_type idx = 0; idx != gsize; ++idx) {
          dFbar_dF0(2, idx) = 0;
          dFbar_dF0(idx, 2) = 0;
          dFbar_dF(2, idx) = 0;
          dFbar_dF(idx, 2) = 0;
        }
        dFbar_dF0(2, 2) = 1;
        dFbar_dF(2, 2) = 1;
      }
      // tangent operator
      const auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
      const auto Kb =
          tfel::math::map<tfel::math::t2tot2<dimension, real>>(Kip.data());
      const auto K_F0 = tfel::math::eval(Kb * dFbar_dF0);
      const auto K_F = tfel::math::eval(Kb * dFbar_dF);
      const auto K_F0v = std::span<const real>(K_F0.data(), gsize * gsize);
      const auto K_Fv = std::span<const real>(K_F.data(), gsize * gsize);
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        if constexpr (evaluateShapeFunctions) {
          child.updateStiffnessMatrix(Ke, K_F0v, shape0, dshape, dshape0, w,
                                      ni);
          child.updateStiffnessMatrix(Ke, K_Fv, shape, dshape, dshape, w, ni);
        } else {
          child.updateStiffnessMatrix(Ke, K_F0v, dshape, dshape0, w, ni);
          child.updateStiffnessMatrix(Ke, K_Fv, dshape, dshape, w, ni);
        }
      }
    }
  }  // end of implementUpdateJacobian

  template <typename Child, Hypothesis H>
  FBarBehaviourIntegratorCRTPBase<Child,
                                  H>::~FBarBehaviourIntegratorCRTPBase() =
      default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FBARBEHAVIOURINTEGRATORCRTPBASE_IXX */
