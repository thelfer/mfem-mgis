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

  template <typename Child>
  void SmallStrainMechanicalBehaviourIntegratorCRTPBase<
      Child>::implementComputeInnerForces(mfem::Vector &Fe,
                                          const mfem::FiniteElement &e,
                                          mfem::ElementTransformation &tr,
                                          const mfem::Vector &u) {
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    // element offset
    const auto nnodes = e.GetDof();
    const auto gsize = this->s1.gradients_stride;
    const auto thsize = this->s1.thermodynamic_forces_stride;
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe.SetSize(e.GetDof() * e.GetDim());
    Fe = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i < ir.GetNPoints(); ++i) {
      // get the gradients of the shape functions
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = ip.weight * tr.Weight();
      // offset of the integration point
      const auto o = eoffset + i;
      auto g = this->s1.gradients.subspan(o * gsize, gsize);
      std::copy(this->macroscopic_gradients.begin(),
                this->macroscopic_gradients.end(), g.begin());
      for (size_type ni = 0; ni != nnodes; ++ni) {
        static_cast<Child *>(this)->updateGradients(g, u, dshape, ni);
      }
      this->integrate(eoffset + i);
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      // assembly of the inner forces
      for (size_type ni = 0; ni != nnodes; ++ni) {
        static_cast<const Child *>(this)->updateInnerForces(Fe, s, dshape, w,
                                                            ni);
      }
    }
  }  // end of implementComputeInnerForces

  template <typename Child>
  void SmallStrainMechanicalBehaviourIntegratorCRTPBase<
      Child>::implementComputeStiffnessMatrix(mfem::DenseMatrix &Ke,
                                              const mfem::FiniteElement &e,
                                              mfem::ElementTransformation &tr) {
#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(e.GetDof(), e.GetDim());
#else
    dshape.SetSize(e.GetDof(), e.GetDim());
#endif
    // element offset
    const auto nnodes = e.GetDof();
    const auto gsize = this->s1.gradients_stride;
    const auto thsize = this->s1.thermodynamic_forces_stride;
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Ke.SetSize(e.GetDof() * e.GetDim(),e.GetDof() * e.GetDim());
    Ke = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i < ir.GetNPoints(); ++i) {
      // get the gradients of the shape functions
      const auto &ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = ip.weight * tr.Weight();
      // offset of the integration point
      const auto o = eoffset + i;
      const auto Kip = this->K.subspan(o * gsize * thsize, gsize * thsize);
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        static_cast<const Child *>(this)->updateStiffnessMatrix(Ke, Kip, dshape,
                                                                w, ni);
      }
    }
  }  // end of implementComputeStiffnessMatrix

  template <typename Child>
  SmallStrainMechanicalBehaviourIntegratorCRTPBase<
      Child>::~SmallStrainMechanicalBehaviourIntegratorCRTPBase() = default;

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<H>::
      SmallStrainMechanicalBehaviourIntegrator(
          const mfem::FiniteElementSpace &fs,
          const size_type m,
          std::shared_ptr<const Behaviour> b_ptr)
      : SmallStrainMechanicalBehaviourIntegratorCRTPBase<
            SmallStrainMechanicalBehaviourIntegrator>(fs, m, std::move(b_ptr)) {
    this->checkHypotheses(H);
  }  // end of SmallStrainMechanicalBehaviourIntegrator

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<
      H>::~SmallStrainMechanicalBehaviourIntegrator() = default;

  /* Specalisations of the methods of the
   * SmallStrainMechanicalBehaviourIntegrator class */

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

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateGradients(
      mgis::span<real>& g,
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
      const auto ux = u[ni];
      const auto uy = u[ni + nnodes];
      const auto uz = u[ni + 2 * nnodes];
      g[0] += ux * dN(ni, 0);                             // xx
      g[1] += uy * dN(ni, 1);                             // yy
      g[2] += uz * dN(ni, 2);                             // zz
      g[3] += (ux * dN(ni, 1) + uy * dN(ni, 0)) * icste;  // xy
      g[4] += (ux * dN(ni, 2) + uz * dN(ni, 0)) * icste;  // xz
      g[5] += (uy * dN(ni, 2) + uz * dN(ni, 1)) * icste;  // yz
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
      const auto nx = ni;
      const auto ny = ni + nnodes;
      const auto nz = ni + 2 * nnodes;
      real B[6][3] = {{dN(ni, 0), 0, 0},                          //
                      {0, dN(ni, 1), 0},                          //
                      {0, 0, dN(ni, 2)},                          //
                      {dN(ni, 1) * icste, dN(ni, 0) * icste, 0},  //
                      {dN(ni, 2) * icste, 0, dN(ni, 0) * icste},  //
                      {0, dN(ni, 2) * icste, dN(ni, 0) * icste}};
      Fe[nx] += w * (B[0][0] * s[0] + B[1][0] * s[1] + B[2][0] * s[2] +
                     B[3][0] * s[3] + B[4][0] * s[4] + B[5][0] * s[5]);
      Fe[ny] += w * (B[0][1] * s[0] + B[1][1] * s[1] + B[2][1] * s[2] +
                     B[3][1] * s[3] + B[4][1] * s[4] + B[5][1] * s[5]);
      Fe[nz] += w * (B[0][2] * s[0] + B[1][2] * s[1] + B[2][2] * s[2] +
                     B[3][2] * s[3] + B[4][2] * s[4] + B[5][2] * s[5]);
    }
  }  // end of updateInnerForces

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateStiffnessMatrix(
      mfem::DenseMatrix &Ke,
      const mgis::span<const real>& Kip,
      const mfem::DenseMatrix &dN,
      const real w,
      const size_type ni) const {
    constexpr const auto icste = 0.70710678118654752440;
    static_assert(H == Hypothesis::TRIDIMENSIONAL, "unsupported hypothesis");
    if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = dN.NumRows();
      real Bi[6][3] = {{dN(ni, 0), 0, 0},                          //
                       {0, dN(ni, 1), 0},                          //
                       {0, 0, dN(ni, 2)},                          //
                       {dN(ni, 1) * icste, dN(ni, 0) * icste, 0},  //
                       {dN(ni, 2) * icste, 0, dN(ni, 0) * icste},  //
                       {0, dN(ni, 2) * icste, dN(ni, 0) * icste}};
      for (size_type nj = 0; nj != nnodes; ++nj) {
        real Bj[6][3] = {{dN(nj, 0), 0, 0},                          //
                         {0, dN(nj, 1), 0},                          //
                         {0, 0, dN(nj, 2)},                          //
                         {dN(nj, 1) * icste, dN(nj, 0) * icste, 0},  //
                         {dN(nj, 2) * icste, 0, dN(nj, 0) * icste},  //
                         {0, dN(nj, 2) * icste, dN(nj, 0) * icste}};
        real KB[6][3];
        for (size_type i = 0; i != 6; ++i) {
          for (size_type j = 0; j != 3; ++j) {
            KB[i][j] = real{};
            for (size_type k = 0; k != 6; ++k) {
              KB[i][j] += Kip[i * 6 + k] * Bj[k][j];
            }
          }
        }
        for (size_type i = 0; i != 3; ++i) {
          for (size_type j = 0; j != 3; ++j) {
            auto tBKB = real{};
            for (size_type k = 0; k != 6; ++k) {
              tBKB += Bi[k][i] * KB[k][j];
            }
            Ke(ni + i * nnodes, nj + j * nnodes) += w * tBKB;
          }
        }
      } // end of for (size_type nj = 0; nj != nnodes; ++nj)
    } // end of if constexpr (H == Hypothesis::TRIDIMENSIONAL)
  }   // end of updateStiffnessMatrix

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX */
