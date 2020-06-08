/*!
 * \file   src/MGISIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <algorithm>
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fem.hpp"
#include "mfem/fem/nonlininteg.hpp"
#include "mfem/fem/fespace.hpp"
#include "MFEMMGIS/MGISIntegrator.hxx"

namespace mfem_mgis {

  void MGISIntegrator::AssembleElementMatrix(const mfem::FiniteElement &el,
                                             mfem::ElementTransformation &Trans,
                                             mfem::DenseMatrix &elmat) {
    int dof = el.GetDof();
    int dim = el.GetDim();
    double w, L, M;

    MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
    Vector divshape(dim * dof);
#else
    dshape.SetSize(dof, dim);
    gshape.SetSize(dof, dim);
    pelmat.SetSize(dof);
    divshape.SetSize(dim * dof);
#endif

    elmat.SetSize(dof * dim);


    const auto *ir = this->IntRule;
    if (ir == NULL) {
      int order = 2 * Trans.OrderGrad(&el);  // correct order?
      ir = &mfem::IntRules.Get(el.GetGeomType(), order);
    }

    elmat = 0.0;

    for (int pt = 0; pt < ir->GetNPoints(); ++pt) {
      const auto &ip = ir->IntPoint(pt);

      el.CalcDShape(ip, dshape);

      Trans.SetIntPoint(&ip);
      w = ip.weight * Trans.Weight();
      Mult(dshape, Trans.InverseJacobian(), gshape);
      MultAAt(gshape, pelmat);
      gshape.GradToDiv(divshape);

      M = mu->Eval(Trans, ip);
      if (lambda) {
        L = lambda->Eval(Trans, ip);
      } else {
        L = q_lambda * M;
        M = q_mu * M;
      }

      if (L != 0.0) {
        AddMult_a_VVt(L * w, divshape, elmat);
      }

      if (M != 0.0) {
        for (int d = 0; d < dim; d++) {
          for (int k = 0; k < dof; k++)
            for (int l = 0; l < dof; l++) {
              elmat(dof * d + k, dof * d + l) += (M * w) * pelmat(k, l);
            }
        }
        for (int i = 0; i < dim; i++)
          for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dof; k++)
              for (int l = 0; l < dof; l++) {
                elmat(dof * i + k, dof * j + l) += (M * w) * gshape(k, j) * gshape(l, i);
              }
          }
      }
    }
  }

  void MGISIntegrator::ComputeElementFlux(const mfem::FiniteElement &el,
                                          mfem::ElementTransformation &Trans,
                                          mfem::Vector &u,
                                          const mfem::FiniteElement &fluxelem,
                                          mfem::Vector &flux,
                                          bool) {
    const int dof = el.GetDof();
    const int dim = el.GetDim();
    const int tdim = dim * (dim + 1) / 2;  // num. entries in a symmetric tensor
    double L, M;

    MFEM_ASSERT(dim == 2 || dim == 3, "dimension is not supported: dim = " << dim);
    MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");
    MFEM_ASSERT(fluxelem.GetMapType() == FiniteElement::VALUE, "");
    MFEM_ASSERT(dynamic_cast<const NodalFiniteElement *>(&fluxelem), "");

#ifdef MFEM_THREAD_SAFE
    DenseMatrix dshape(dof, dim);
#else
    dshape.SetSize(dof, dim);
#endif

    double gh_data[9], grad_data[9];
    mfem::DenseMatrix gh(gh_data, dim, dim);
    mfem::DenseMatrix grad(grad_data, dim, dim);

    const auto&ir = fluxelem.GetNodes();
    const int fnd = ir.GetNPoints();
    flux.SetSize(fnd * tdim);

    mfem::DenseMatrix loc_data_mat(u.GetData(), dof, dim);
    for (int i = 0; i < fnd; i++) {
      const auto&ip = ir.IntPoint(i);
      el.CalcDShape(ip, dshape);
      MultAtB(loc_data_mat, dshape, gh);

      Trans.SetIntPoint(&ip);
      Mult(gh, Trans.InverseJacobian(), grad);

      M = mu->Eval(Trans, ip);
      if (lambda) {
        L = lambda->Eval(Trans, ip);
      } else {
        L = q_lambda * M;
        M = q_mu * M;
      }

      // stress = 2*M*e(u) + L*tr(e(u))*I, where
      //   e(u) = (1/2)*(grad(u) + grad(u)^T)
      const double M2 = 2.0 * M;
      if (dim == 2) {
        L *= (grad(0, 0) + grad(1, 1));
        // order of the stress entries: s_xx, s_yy, s_xy
        flux(i + fnd * 0) = M2 * grad(0, 0) + L;
        flux(i + fnd * 1) = M2 * grad(1, 1) + L;
        flux(i + fnd * 2) = M * (grad(0, 1) + grad(1, 0));
      } else if (dim == 3) {
        L *= (grad(0, 0) + grad(1, 1) + grad(2, 2));
        // order of the stress entries: s_xx, s_yy, s_zz, s_xy, s_xz, s_yz
        flux(i + fnd * 0) = M2 * grad(0, 0) + L;
        flux(i + fnd * 1) = M2 * grad(1, 1) + L;
        flux(i + fnd * 2) = M2 * grad(2, 2) + L;
        flux(i + fnd * 3) = M * (grad(0, 1) + grad(1, 0));
        flux(i + fnd * 4) = M * (grad(0, 2) + grad(2, 0));
        flux(i + fnd * 5) = M * (grad(1, 2) + grad(2, 1));
      }
    }
  }

  double MGISIntegrator::ComputeFluxEnergy(const mfem::FiniteElement &fluxelem,
                                           mfem::ElementTransformation &Trans,
                                           mfem::Vector &flux,
                                           mfem::Vector *) {
    const int dof = fluxelem.GetDof();
    const int dim = fluxelem.GetDim();
    const int tdim = dim * (dim + 1) / 2;  // num. entries in a symmetric tensor
    double L, M;

    // The MFEM_ASSERT constraints in MGISIntegrator::ComputeElementFlux
    // are assumed here too.
    MFEM_ASSERT(d_energy == NULL, "anisotropic estimates are not supported");
    MFEM_ASSERT(flux.Size() == dof * tdim, "invalid 'flux' vector");

#ifndef MFEM_THREAD_SAFE
    shape.SetSize(dof);
#else
    Vector shape(dof);
#endif
    double pointstress_data[6];
    mfem::Vector pointstress(pointstress_data, tdim);

    // View of the 'flux' vector as a (dof x tdim) matrix
    mfem::DenseMatrix flux_mat(flux.GetData(), dof, tdim);

    // Use the same integration rule as in AssembleElementMatrix, replacing 'el'
    // with 'fluxelem' when 'IntRule' is not set.
    // Should we be using a different (more accurate) rule here?
    const auto*ir = this->IntRule;
    if (ir == NULL) {
      int order = 2 * Trans.OrderGrad(&fluxelem);
      ir = &mfem::IntRules.Get(fluxelem.GetGeomType(), order);
    }

    double energy = 0.0;

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const auto &ip = ir->IntPoint(i);
      fluxelem.CalcShape(ip, shape);

      flux_mat.MultTranspose(shape, pointstress);

      Trans.SetIntPoint(&ip);
      double w = Trans.Weight() * ip.weight;

      M = mu->Eval(Trans, ip);
      if (lambda) {
        L = lambda->Eval(Trans, ip);
      } else {
        L = q_lambda * M;
        M = q_mu * M;
      }

      // The strain energy density at a point is given by (1/2)*(s : e) where s
      // and e are the stress and strain tensors, respectively. Since we only
      // have the stress, we need to compute the strain from the stress:
      //    s = 2*mu*e + lambda*tr(e)*I
      // Taking trace on both sides we find:
      //    tr(s) = 2*mu*tr(e) + lambda*tr(e)*dim = (2*mu + dim*lambda)*tr(e)
      // which gives:
      //    tr(e) = tr(s)/(2*mu + dim*lambda)
      // Then from the first identity above we can find the strain:
      //    e = (1/(2*mu))*(s - lambda*tr(e)*I)

      double pt_e;  // point strain energy density
      const double *s = pointstress_data;
      if (dim == 2) {
        // s entries: s_xx, s_yy, s_xy
        const double tr_e = (s[0] + s[1]) / (2 * (M + L));
        L *= tr_e;
        pt_e = (0.25 / M) * (s[0] * (s[0] - L) + s[1] * (s[1] - L) + 2 * s[2] * s[2]);
      } else  // (dim == 3)
      {
        // s entries: s_xx, s_yy, s_zz, s_xy, s_xz, s_yz
        const double tr_e = (s[0] + s[1] + s[2]) / (2 * M + 3 * L);
        L *= tr_e;
        pt_e = (0.25 / M) * (s[0] * (s[0] - L) + s[1] * (s[1] - L) + s[2] * (s[2] - L) +
                             2 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
      }

      energy += w * pt_e;
    }

    return energy;
  }

  MGISIntegrator::~MGISIntegrator() = default;

}  // end of namespace mfem_mgis
