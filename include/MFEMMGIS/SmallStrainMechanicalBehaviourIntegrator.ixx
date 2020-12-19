/*!
 * \file   SmallStrainMechanicalBehaviourIntegrator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#ifndef LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX
#define LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX

#include <utility>
#include "mfem/fem/eltrans.hpp"
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  template <Hypothesis H>
  const mfem::IntegrationRule &
  SmallStrainMechanicalBehaviourIntegrator<H>::getIntegrationRule(
      const mfem::FiniteElement &el, const mfem::ElementTransformation &Trans) {
    const auto order = 2 * Trans.OrderGrad(&el);
    if constexpr ((H == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) ||
                  (H == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS) ||
                  (H == Hypothesis::AXISYMMETRICAL)) {
      const auto &ir = mfem::IntRules.Get(el.GetGeomType(), order + 1);
      return ir;
    } else {
      const auto &ir = mfem::IntRules.Get(el.GetGeomType(), order);
      return ir;
    }
  }  // end of getIntegrationRule

  template <Hypothesis H>
  std::shared_ptr<const PartialQuadratureSpace>
  SmallStrainMechanicalBehaviourIntegrator<H>::buildQuadratureSpace(
      const FiniteElementDiscretization &fed, const size_type m) {
    auto selector = [](const mfem::FiniteElement &e,
                       const mfem::ElementTransformation &tr)
        -> const mfem::IntegrationRule & {
      return getIntegrationRule(e, tr);
    };  // end of selector
    return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
  }  // end of buildQuadratureSpace

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<H>::
      SmallStrainMechanicalBehaviourIntegrator(
          const FiniteElementDiscretization &fed,
          const size_type m,
          std::unique_ptr<const Behaviour> b_ptr)
      : StandardBehaviourIntegratorCRTPBase<
            SmallStrainMechanicalBehaviourIntegrator>(
            buildQuadratureSpace(fed, m), std::move(b_ptr)) {
    this->checkHypotheses(H);
    if (this->b.btype != Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      this->throwInvalidBehaviourType(
          "SmallStrainMechanicalBehaviourIntegrator",
          "expected a strain based behaviour");
    }
    if (this->b.kinematic != Behaviour::SMALLSTRAINKINEMATIC) {
      this->throwInvalidBehaviourKinematic(
          "SmallStrainMechanicalBehaviourIntegrator",
          "expected a small strain based behaviour");
    }
  }  // end of SmallStrainMechanicalBehaviourIntegrator

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<
      H>::~SmallStrainMechanicalBehaviourIntegrator() = default;

  // inline implementations

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateGradients(
      mgis::span<real> &g,
      const mfem::Vector &u,
      const mfem::DenseMatrix &dN,
      const size_type ni) {
    constexpr const auto icste = 0.70710678118654752440;
    static_assert((H == Hypothesis::TRIDIMENSIONAL) ||
                      (H == Hypothesis::PLANESTRAIN) ||
                      (H == Hypothesis::PLANESTRESS),
                  "unsupported hypothesis");
    if constexpr ((H == Hypothesis::PLANESTRAIN) ||
                  (H == Hypothesis::PLANESTRESS)) {
      const auto nnodes = dN.NumRows();
      const auto nx = ni;
      const auto ny = ni + nnodes;
      const auto ux = u[nx];
      const auto uy = u[ny];
      g[0] += ux * dN(nx, 0);                             // xx
      g[1] += uy * dN(ny, 1);                             // yy
      g[2] += 0;                                          // zz
      g[3] += (ux * dN(nx, 1) + uy * dN(ny, 0)) * icste;  // xy
    } else if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = dN.NumRows();
      const auto Bi_0_0 = dN(ni, 0);
      const auto Bi_1_1 = dN(ni, 1);
      const auto Bi_2_2 = dN(ni, 2);
      const auto Bi_3_0 = dN(ni, 1) * icste;
      const auto Bi_3_1 = dN(ni, 0) * icste;
      const auto Bi_4_0 = dN(ni, 2) * icste;
      const auto Bi_4_2 = dN(ni, 0) * icste;
      const auto Bi_5_1 = dN(ni, 2) * icste;
      const auto Bi_5_2 = dN(ni, 1) * icste;
      const auto u_0 = u[ni];
      const auto u_1 = u[ni + nnodes];
      const auto u_2 = u[ni + 2 * nnodes];
      g[0] += u_0 * Bi_0_0;
      g[1] += u_1 * Bi_1_1;
      g[2] += u_2 * Bi_2_2;
      g[3] += u_0 * Bi_3_0 + u_1 * Bi_3_1;
      g[4] += u_2 * Bi_4_2 + Bi_4_0 * u_0;
      g[5] += u_2 * Bi_5_2 + u_1 * Bi_5_1;
    }
  }  // end of updateGradients

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateInnerForces(
      mfem::Vector &Fe,
      const mgis::span<const real> &s,
      const mfem::DenseMatrix &dN,
      const real w,
      const size_type ni) const {
    constexpr const auto icste = 0.70710678118654752440;
    static_assert(H == Hypothesis::TRIDIMENSIONAL, "unsupported hypothesis");
    if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = dN.NumRows();
      const auto Bi_0_0 = dN(ni, 0);
      const auto Bi_1_1 = dN(ni, 1);
      const auto Bi_2_2 = dN(ni, 2);
      const auto Bi_3_0 = dN(ni, 1) * icste;
      const auto Bi_3_1 = dN(ni, 0) * icste;
      const auto Bi_4_0 = dN(ni, 2) * icste;
      const auto Bi_4_2 = dN(ni, 0) * icste;
      const auto Bi_5_1 = dN(ni, 2) * icste;
      const auto Bi_5_2 = dN(ni, 1) * icste;
      const auto ni_0 = ni;
      const auto ni_1 = ni + nnodes;
      const auto ni_2 = ni + 2 * nnodes;
      Fe[ni_0] += w * s[4] * Bi_4_0 + Bi_3_0 * w * s[3] + w * Bi_0_0 * s[0];
      Fe[ni_1] += w * Bi_5_1 * s[5] + w * Bi_3_1 * s[3] + w * s[1] * Bi_1_1;
      Fe[ni_2] += w * s[5] * Bi_5_2 + w * s[2] * Bi_2_2 + w * s[4] * Bi_4_2;
    }
  }  // end of updateInnerForces

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateStiffnessMatrix(
      mfem::DenseMatrix &Ke,
      const mgis::span<const real> &Kip,
      const mfem::DenseMatrix &dN,
      const real w,
      const size_type ni) const {
    constexpr const auto icste = 0.70710678118654752440;
    static_assert(H == Hypothesis::TRIDIMENSIONAL, "unsupported hypothesis");
    if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = dN.NumRows();
      const auto Bi_0_0 = dN(ni, 0);
      const auto Bi_1_1 = dN(ni, 1);
      const auto Bi_2_2 = dN(ni, 2);
      const auto Bi_3_0 = dN(ni, 1) * icste;
      const auto Bi_3_1 = dN(ni, 0) * icste;
      const auto Bi_4_0 = dN(ni, 2) * icste;
      const auto Bi_4_2 = dN(ni, 0) * icste;
      const auto Bi_5_1 = dN(ni, 2) * icste;
      const auto Bi_5_2 = dN(ni, 1) * icste;
      const auto ni_0 = ni;
      const auto ni_1 = ni + nnodes;
      const auto ni_2 = ni + 2 * nnodes;
      for (size_type nj = 0; nj != nnodes; ++nj) {
        const auto Bj_0_0 = dN(nj, 0);
        const auto Bj_1_1 = dN(nj, 1);
        const auto Bj_2_2 = dN(nj, 2);
        const auto Bj_3_0 = dN(nj, 1) * icste;
        const auto Bj_3_1 = dN(nj, 0) * icste;
        const auto Bj_4_0 = dN(nj, 2) * icste;
        const auto Bj_4_2 = dN(nj, 0) * icste;
        const auto Bj_5_1 = dN(nj, 2) * icste;
        const auto Bj_5_2 = dN(nj, 1) * icste;
        const auto nj_0 = nj;
        const auto nj_1 = nj + nnodes;
        const auto nj_2 = nj + 2 * nnodes;
        Ke(ni_0, nj_0) +=
            Kip[27] * Bi_4_0 * Bj_3_0 * w + Kip[21] * Bi_3_0 * Bj_3_0 * w +
            Kip[18] * Bi_3_0 * w * Bj_0_0 + Kip[3] * Bi_0_0 * Bj_3_0 * w +
            Kip[24] * Bi_4_0 * w * Bj_0_0 + Bi_0_0 * Kip[0] * w * Bj_0_0 +
            Bi_4_0 * Kip[28] * Bj_4_0 * w + Kip[22] * Bj_4_0 * Bi_3_0 * w +
            Bj_4_0 * Bi_0_0 * Kip[4] * w;
        Ke(ni_0, nj_1) +=
            Bj_3_1 * Kip[3] * Bi_0_0 * w + Bj_5_1 * Kip[5] * Bi_0_0 * w +
            Bj_5_1 * Bi_4_0 * Kip[29] * w + Bj_1_1 * Bi_0_0 * Kip[1] * w +
            Kip[21] * Bj_3_1 * Bi_3_0 * w + Bj_1_1 * Bi_4_0 * Kip[25] * w +
            Bj_5_1 * Kip[23] * Bi_3_0 * w + Kip[27] * Bj_3_1 * Bi_4_0 * w +
            Bj_1_1 * Kip[19] * Bi_3_0 * w;
        Ke(ni_0, nj_2) +=
            Bi_4_0 * Bj_5_2 * Kip[29] * w + Bj_2_2 * Bi_3_0 * Kip[20] * w +
            Kip[2] * Bj_2_2 * Bi_0_0 * w + Bi_0_0 * Bj_4_2 * Kip[4] * w +
            Bj_5_2 * Kip[23] * Bi_3_0 * w + Bi_4_0 * Bj_2_2 * Kip[26] * w +
            Kip[22] * Bi_3_0 * Bj_4_2 * w + Kip[5] * Bj_5_2 * Bi_0_0 * w +
            Bi_4_0 * Kip[28] * Bj_4_2 * w;
        Ke(ni_1, nj_0) +=
            Kip[34] * Bj_4_0 * w * Bi_5_1 + Bi_1_1 * Kip[6] * w * Bj_0_0 +
            Kip[33] * Bj_3_0 * w * Bi_5_1 + Bi_3_1 * Kip[22] * Bj_4_0 * w +
            Bi_3_1 * Kip[21] * Bj_3_0 * w + Bi_1_1 * Kip[9] * Bj_3_0 * w +
            Bi_1_1 * Bj_4_0 * Kip[10] * w + Bi_3_1 * Kip[18] * w * Bj_0_0 +
            Kip[30] * w * Bj_0_0 * Bi_5_1;
        Ke(ni_1, nj_1) +=
            Bi_1_1 * Bj_3_1 * Kip[9] * w + Bi_1_1 * Bj_1_1 * Kip[7] * w +
            Bi_3_1 * Kip[21] * Bj_3_1 * w + Kip[11] * Bj_5_1 * Bi_1_1 * w +
            Kip[33] * Bj_3_1 * w * Bi_5_1 + Bi_3_1 * Bj_5_1 * Kip[23] * w +
            Bj_1_1 * Kip[31] * w * Bi_5_1 + Bj_5_1 * Kip[35] * w * Bi_5_1 +
            Bi_3_1 * Bj_1_1 * Kip[19] * w;
        Ke(ni_1, nj_2) +=
            Bj_2_2 * Kip[32] * w * Bi_5_1 + Kip[8] * Bi_1_1 * Bj_2_2 * w +
            Bj_5_2 * Kip[35] * w * Bi_5_1 + Bi_3_1 * Kip[22] * Bj_4_2 * w +
            Kip[11] * Bi_1_1 * Bj_5_2 * w + Kip[34] * Bj_4_2 * w * Bi_5_1 +
            Bi_3_1 * Bj_2_2 * Kip[20] * w + Bi_3_1 * Bj_5_2 * Kip[23] * w +
            Bi_1_1 * Kip[10] * Bj_4_2 * w;
        Ke(ni_2, nj_0) +=
            Kip[24] * Bi_4_2 * w * Bj_0_0 + Kip[12] * Bi_2_2 * w * Bj_0_0 +
            Kip[33] * Bi_5_2 * Bj_3_0 * w + Bi_2_2 * Bj_4_0 * Kip[16] * w +
            Kip[30] * Bi_5_2 * w * Bj_0_0 + Kip[28] * Bj_4_0 * Bi_4_2 * w +
            Kip[15] * Bi_2_2 * Bj_3_0 * w + Bi_5_2 * Kip[34] * Bj_4_0 * w +
            Kip[27] * Bi_4_2 * Bj_3_0 * w;
        Ke(ni_2, nj_1) +=
            Bj_1_1 * Bi_2_2 * Kip[13] * w + Kip[33] * Bj_3_1 * Bi_5_2 * w +
            Bj_1_1 * Kip[25] * Bi_4_2 * w + Bj_5_1 * Bi_2_2 * Kip[17] * w +
            Bj_5_1 * Bi_5_2 * Kip[35] * w + Bj_3_1 * Kip[15] * Bi_2_2 * w +
            Bj_5_1 * Bi_4_2 * Kip[29] * w + Bj_1_1 * Bi_5_2 * Kip[31] * w +
            Kip[27] * Bj_3_1 * Bi_4_2 * w;
        Ke(ni_2, nj_2) +=
            Bj_5_2 * Bi_4_2 * Kip[29] * w + Bi_2_2 * Kip[16] * Bj_4_2 * w +
            Bi_5_2 * Kip[34] * Bj_4_2 * w + Bi_5_2 * Bj_5_2 * Kip[35] * w +
            Bj_5_2 * Bi_2_2 * Kip[17] * w + Kip[28] * Bi_4_2 * Bj_4_2 * w +
            Bi_2_2 * Bj_2_2 * Kip[14] * w + Bi_5_2 * Bj_2_2 * Kip[32] * w +
            Bj_2_2 * Bi_4_2 * Kip[26] * w;
        }  // end of for (size_type nj = 0; nj != nnodes; ++nj)
      }  // end of if constexpr (H == Hypothesis::TRIDIMENSIONAL)
    }    // end of updateStiffnessMatrix

  /* Methods that must be explicitely instanciated in a source file */

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeInnerForces(mfem::Vector &,
                         const mfem::FiniteElement &,
                         mfem::ElementTransformation &,
                         const mfem::Vector &);

  template <>
  void SmallStrainMechanicalBehaviourIntegrator<Hypothesis::TRIDIMENSIONAL>::
      computeStiffnessMatrix(mfem::DenseMatrix &,
                             const mfem::FiniteElement &,
                             mfem::ElementTransformation &,
                             const mfem::Vector &);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX */
