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
    const auto nnodes = e.GetDof() / e.GetDim();
    const auto gsize = this->s1.gradients_stride;
    const auto thsize = this->s1.thermodynamic_forces_stride;
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Fe = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i < ir.GetNPoints(); i++) {
      // get the gradients of the shape functions
      const auto&ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape);
	// get the weights associated to point ip
      const auto w = ip.weight * tr.Weight();
      // offset of the integration point
      const auto o = eoffset + i;
      auto g = this->s1.gradients.subspan(o * gsize, gsize);
      std::fill(g.begin(), g.end(), real(0));
      for (size_type ni = 0; ni != nnodes; ++ni) {
        static_cast<Child *>(this)->updateStrain(g.data(), u, dshape, ni);
      }
      this->integrate(eoffset + i);
      const auto s = this->s1.thermodynamic_forces.subspan(o * thsize, thsize);
      // assembly of the inner forces
      for (size_type ni = 0; ni != nnodes; ++ni) {
        static_cast<const Child *>(this)->updateInnerForces(Fe, s.data(),
                                                            dshape, w, ni);
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
    const auto nnodes = e.GetDof() / e.GetDim();
    const auto gsize = this->s1.gradients_stride;
    const auto thsize = this->s1.thermodynamic_forces_stride;
    const auto eoffset = this->quadrature_space->getOffset(tr.ElementNo);
    Ke = 0.;
    const auto &ir = this->getIntegrationRule(e, tr);
    for (size_type i = 0; i < ir.GetNPoints(); i++) {
      // get the gradients of the shape functions
      const auto&ip = ir.IntPoint(i);
      tr.SetIntPoint(&ip);
      e.CalcPhysDShape(tr, dshape);
      // get the weights associated to point ip
      const auto w = ip.weight * tr.Weight();
      // offset of the integration point
      const auto o = eoffset + i;
      auto Kip = this->K.subspan(o * gsize * thsize, gsize * thsize);
      std::fill(Kip.begin(), Kip.end(), real(0));
      // assembly of the stiffness matrix
      for (size_type ni = 0; ni != nnodes; ++ni) {
        for (size_type nj = 0; nj != nnodes; ++nj) {

        }
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
          std::unique_ptr<const Behaviour> b_ptr)
      : SmallStrainMechanicalBehaviourIntegratorCRTPBase<
            SmallStrainMechanicalBehaviourIntegrator>(fs, m, std::move(b_ptr)) {
    this->checkHypotheses(H);
  }  // end of SmallStrainMechanicalBehaviourIntegrator

  template <Hypothesis H>
  SmallStrainMechanicalBehaviourIntegrator<
      H>::~SmallStrainMechanicalBehaviourIntegrator() = default;

  /* Specalisations of the methods of the  SmallStrainMechanicalBehaviourIntegrator class */

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
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateStrain(
      real *const g,
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
      const auto nnodes = u.Size() / 2;
      const auto nx = ni;
      const auto ny = ni + nnodes;
      const auto ux = u[nx];
      const auto uy = u[ny];
      g[0] += ux * dN(nx, 0);                             // xx
      g[1] += uy * dN(ny, 1);                             // yy
      g[2] += 0;                                          // zz
      g[3] += (ux * dN(nx, 1) + uy * dN(ny, 0)) * icste;  // xy
    } else if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = u.Size() / 3;
      const auto nx = ni;
      const auto ny = ni + nnodes;
      const auto nz = ni + 2 * nnodes;
      const auto ux = u[nx];
      const auto uy = u[ny];
      const auto uz = u[nz];
      g[0] += ux * dN(nx, 0);                             // xx
      g[1] += uy * dN(ny, 1);                             // yy
      g[2] += uz * dN(nz, 2);                             // zz
      g[3] += (ux * dN(nx, 1) + uy * dN(ny, 0)) * icste;  // xy
      g[4] += (ux * dN(nx, 2) + uz * dN(nz, 0)) * icste;  // xz
      g[5] += (uy * dN(ny, 2) + uz * dN(nz, 1)) * icste;  // yz
    }
  }  // end of updateStrain

  template <Hypothesis H>
  void SmallStrainMechanicalBehaviourIntegrator<H>::updateInnerForces(
      mfem::Vector &Fe,
      const real *const s,
      const mfem::DenseMatrix &dN,
      const real w,
      const size_type ni) const {
    constexpr const auto icste = 0.70710678118654752440;
    static_assert((H == Hypothesis::TRIDIMENSIONAL) ||
                      (H == Hypothesis::PLANESTRAIN) ||
                      (H == Hypothesis::PLANESTRESS),
                  "unsupported hypothesis");
    if constexpr (H == Hypothesis::TRIDIMENSIONAL) {
      const auto nnodes = dN.NumRows() / 3;
      const auto nx = ni;
      const auto ny = ni + nnodes;
      const auto nz = ni + 2 * nnodes;
      real B[6][3] = {{dN(nx, 0), 0, 0},                          //
                      {0, dN(ny, 0), 0},                          //
                      {0, 0, dN(nz, 0)},                          //
                      {dN(nx, 1) * icste, dN(ny, 0) * icste, 0},  //
                      {dN(nx, 2) * icste, 0, dN(nz, 0) * icste},  //
                      {0, dN(ny, 2) * icste, dN(nz, 0) * icste}};
      Fe[nx] += w * (B[0][0] * s[0] + B[1][0] * s[1] + B[2][0] * s[2] +
                     B[3][0] * s[3] + B[4][0] * s[4] + B[5][0] * s[5]);
      Fe[ny] += w * (B[0][1] * s[0] + B[1][1] * s[1] + B[2][1] * s[2] +
                     B[3][1] * s[3] + B[4][1] * s[4] + B[5][1] * s[5]);
      Fe[nz] += w * (B[0][2] * s[0] + B[1][2] * s[1] + B[2][2] * s[2] +
                     B[3][2] * s[3] + B[4][2] * s[4] + B[5][2] * s[5]);
    }
  }  // end of updateInnerForces

  //     const auto dim = el.GetDim();
  //     const auto dof = el.GetDof();
  //     // size of the gradients
  //     const auto gs = mgis::getStensorSize(H);
  //
  // #ifdef MFEM_MGIS_DEBUG
  //     this->checkDimension();
  //     MFEM_VERIFY(dim == Trans.GetSpaceDim(), "");
  // #endif
  //
  // #ifdef MFEM_THREAD_SAFE
  //     DenseMatrix gshape(dof, dim);
  // #else
  //     Emat.SetSize(gs);
  //     dshape.SetSize(dof, dim);
  //     tBEBmat.SetSize(dim, dim);
  //     tBmatL.SetSize(dim, gs);
  //     tBmatR.SetSize(dim, gs);
  //     BmatR.SetSize(gs, dim);
  //     EBmat.SetSize(gs, dim);
  //     gtBmat.SetSize(dim, dof * gs);
  // #endif
  //
  //     Ke.SetSize(dof * dim);
  //     Ke = 0.0;
  //     //
  //     const auto &ir = this->getIntegrationRule();
  //     for (int i = 0; i < ir -> GetNPoints(); i++){
  //       const auto&ip = ir->IntPoint(i);
  //       tr.SetIntPoint(&ip);
  //       // Each row of the result dshape contains
  //       // the derivatives of one shape function at the point ip.
  //       e.CalcPhysDShape(tr, dshape);
  //       // Get the transformation Trans for point ip
  //       // Get the weights associated to point ip
  //       w = ip.weight * tr.Weight();
  //       // compute the strain at the end of the time step
  //       this->computeStrain(u, dshape);
  //       this->integrate();
  //
  //       // Fill gtBmat
  //       gtBmat = 0.;
  //       for (int s = 0; s < dof; s++) {
  //         int offset = gs * s;
  //         // Fill gtBmat that depends on dof 's'
  //         // gtBmat(:,offset:offset+gs) contains tBmat for dof 's'
  //         for (int k = 0; k < dim; k++) {
  //           gtBmat(k, offset + k) = dshape(s, k);
  //         }
  //         gtBmat(0, offset + dim) = dshape(s, 1);  // N_y
  //         gtBmat(1, offset + dim) = dshape(s, 0);  // N_x
  //         if (dim == 3) {
  //           gtBmat(1, offset + 4) = dshape(s, 2);  // N_z
  //           gtBmat(2, offset + 4) = dshape(s, 1);  // N_y
  //           gtBmat(0, offset + 5) = dshape(s, 2);  // N_z
  //           gtBmat(2, offset + 5) = dshape(s, 0);  // N_x
  //         }
  //       }
  //       // Perform the main matrix multiplications
  //       for (int s = 0; s < dof; s++) {
  //         tBmatL.CopyCols(gtBmat, gs * s, gs * s + gs - 1);
  //         for (int t = 0; t < dof; t++) {
  //           tBmatR.CopyCols(gtBmat, gs * t, gs * t + gs - 1);
  //           BmatR.Transpose(tBmatR);
  //           Mult(Emat, BmatR, EBmat);
  //           Mult(tBmatL, EBmat, tBEBmat);
  //           // Multiply by the weight
  //           tBEBmat *= w;
  //           // Add the contribution to the RHS
  //           for (int k = 0; k < dim; k++){
  //             for (int l = 0; l < dim; l++){
  //               Ke(dof * k + s, dof * l + t) += tBEBmat(k, l);
  //             }
  //           }
  //         }
  //       }
  //     }

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SMALLSTRAINMECHANICALBEHAVIOURINTEGRATOR_IXX */
