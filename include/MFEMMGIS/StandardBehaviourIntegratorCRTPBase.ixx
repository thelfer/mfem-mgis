/*!
 * \file   include/MFEMMGIS/StandardBehaviourIntegratorCRTPBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2020
 */

#ifndef LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX
#define LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX

#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/BehaviourIntegratorTraits.hxx"

namespace mfem_mgis {

  template <typename Child>
  bool StandardBehaviourIntegratorCRTPBase<Child>::implementIntegrate(
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
    } else {
#ifdef MFEM_THREAD_SAFE
      mfem::Vector shape;
      mfem::DenseMatrix dshape;
      if constexpr (evaluateShapeFunctions) {
        shape.SetSize(e.GetDof());
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        dshape.setSize(e.GetDof(), e.GetDim());
      }
#else
      if constexpr (evaluateShapeFunctions) {
        this->shape.SetSize(e.GetDof());
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        this->dshape.SetSize(e.GetDof(), e.GetDim());
      }
#endif
      const auto nnodes = e.GetDof();
      const auto gsize = this->s1.gradients_stride;
      const auto thsize = this->s1.thermodynamic_forces_stride;
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
        if constexpr (evaluateShapeFunctionsDerivatives) {
          // get the gradients of the shape functions
          e.CalcPhysDShape(tr, dshape);
        }
        auto g = this->s1.gradients.subspan(o * gsize, gsize);
        std::copy(this->macroscopic_gradients.begin(),
                  this->macroscopic_gradients.end(), g.begin());
        for (size_type ni = 0; ni != nnodes; ++ni) {
          if constexpr (
              (Traits::gradientsComputationRequiresShapeFunctions) &&
              (Traits::gradientsComputationRequiresShapeFunctionsDerivatives)) {
            child.updateGradients(g, u, shape, dshape, ni);
          } else if constexpr (
              Traits::gradientsComputationRequiresShapeFunctionsDerivatives) {
            child.updateGradients(g, u, dshape, ni);
          } else {
            child.updateGradients(g, u, shape, ni);
          }
        }
        const auto r = child.getRotationMatrix(o);
        child.rotateGradients(g, r);
        if (!this->performsLocalBehaviourIntegration(o, it)) {
          return false;
        }
        const auto s =
            this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
        // Here we rotate the tangent operator blocks but not the thermodynamic
        // forces.
        if (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
          auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
          child.rotateTangentOperatorBlocks(Kip, r);
        }
      }
    }
    return true;
  }  // end of implementIntegrate

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementUpdateResidual(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &) {
    this->implementComputeInnerForces(Fe, e, tr);
  }  // end of implementUpdateResidual

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementComputeInnerForces(
      mfem::Vector &Fe,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
    auto &child = static_cast<Child &>(*this);
#ifdef MFEM_THREAD_SAFE
      mfem::Vector shape;
      mfem::DenseMatrix dshape;
      if constexpr (evaluateShapeFunctions) {
        shape.SetSize(e.GetDof());
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        dshape.setSize(e.GetDof(), e.GetDim());
      }
#else
      if constexpr (evaluateShapeFunctions) {
        this->shape.SetSize(e.GetDof());
      }
      if constexpr (evaluateShapeFunctionsDerivatives) {
        this->dshape.SetSize(e.GetDof(), e.GetDim());
      }
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

  template <typename Child>
  void StandardBehaviourIntegratorCRTPBase<Child>::implementUpdateJacobian(
      mfem::DenseMatrix &Ke,
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr) {
    using Traits = BehaviourIntegratorTraits<Child>;
    auto &child = static_cast<Child &>(*this);
    constexpr const auto evaluateShapeFunctions =
        Traits::updateExternalStateVariablesFromUnknownsValues ||
        Traits::gradientsComputationRequiresShapeFunctions;
    constexpr const auto evaluateShapeFunctionsDerivatives =
        Traits::gradientsComputationRequiresShapeFunctionsDerivatives;
    static_assert(
        ((evaluateShapeFunctions) || (evaluateShapeFunctionsDerivatives)),
        "neither the shape functions or their derivatives are required to "
        "compute the gradients, this is not supported");
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape;
    if constexpr (evaluateShapeFunctions) {
      shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      dshape.SetSize(e.GetDof(), e.GetDim());
    }
  }
#else
    if constexpr (evaluateShapeFunctions) {
      this->shape.SetSize(e.GetDof());
    }
    if constexpr (evaluateShapeFunctionsDerivatives) {
      this->dshape.SetSize(e.GetDof(), e.GetDim());
    }
#endif
  // element offset
  const auto nnodes = e.GetDof();
  const auto eoffset = this -> quadrature_space -> getOffset(tr.ElementNo);
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
    if constexpr (evaluateShapeFunctionsDerivatives) {
      // get the derivatives of the shape functions
      e.CalcPhysDShape(tr, dshape);
    }
    // get the weights associated to point ip
    const auto w = child.getIntegrationPointWeight(tr, ip);
    // offset of the integration point
    const auto o = eoffset + i;
    const auto Kip = this->K.subspan(o * (this->K_stride), this->K_stride);
    // assembly of the stiffness matrix
    for (size_type ni = 0; ni != nnodes; ++ni) {
      if constexpr ((evaluateShapeFunctions) &&
                    (evaluateShapeFunctionsDerivatives)) {
        child.updateStiffnessMatrix(Ke, Kip, shape, dshape, w, ni);
      } else if constexpr (evaluateShapeFunctionsDerivatives) {
        child.updateStiffnessMatrix(Ke, Kip, dshape, w, ni);
      } else {
        child.updateStiffnessMatrix(Ke, Kip, shape, w, ni);
      }
    }
  }
}  // namespace mfem_mgis

template <typename Child>
StandardBehaviourIntegratorCRTPBase<
    Child>::~StandardBehaviourIntegratorCRTPBase() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_STANDARDBEHAVIOURINTEGRATORCRTPBASE_IXX */
