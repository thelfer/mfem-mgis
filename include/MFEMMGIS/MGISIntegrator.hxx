/*!
 * \file   include/MFEMMGIS/MGISIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_MGISINTEGRATOR_HXX
#define LIB_MFEM_MGIS_MGISINTEGRATOR_HXX

#include "mfem/fem/bilinearform.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /** Integrator for the linear elasticity form:
      a(u,v) = (lambda div(u), div(v)) + (2 mu e(u), e(v)),
      where e(v) = (1/2) (grad(v) + grad(v)^T).
      This is a 'Vector' integrator, i.e. defined for FE spaces
      using multiple copies of a scalar FE space. */
  struct MFEM_MGIS_EXPORT MGISIntegrator : public mfem::BilinearFormIntegrator {

    MGISIntegrator(mfem::Coefficient &l, mfem::Coefficient &m) {
      lambda = &l;
      mu = &m;
    }
    /** With this constructor lambda = q_l * m and mu = q_m * m;
        if dim * q_l + 2 * q_m = 0 then trace(sigma) = 0. */
    MGISIntegrator(mfem::Coefficient &m, double q_l, double q_m) {
      lambda = NULL;
      mu = &m;
      q_lambda = q_l;
      q_mu = q_m;
    }

    void AssembleElementMatrix(const mfem::FiniteElement &,
                               mfem::ElementTransformation &,
                               mfem::DenseMatrix &) override;

    /** Compute the stress corresponding to the local displacement @a u and
        interpolate it at the nodes of the given @a fluxelem. Only the symmetric
        part of the stress is stored, so that the size of @a flux is equal to
        the number of DOFs in @a fluxelem times dim*(dim+1)/2. In 2D, the order
        of the stress components is: s_xx, s_yy, s_xy. In 3D, it is: s_xx, s_yy,
        s_zz, s_xy, s_xz, s_yz. In other words, @a flux is the local vector for
        a FE space with dim*(dim+1)/2 vector components, based on the finite
        element @a fluxelem. */
    void ComputeElementFlux(const mfem::FiniteElement &el,
                            mfem::ElementTransformation &Trans,
                            mfem::Vector &u,
                            const mfem::FiniteElement &fluxelem,
                            mfem::Vector &flux,
                            bool with_coef = true) override;

    /** Compute the element energy (integral of the strain energy density)
        corresponding to the stress represented by @a flux which is a vector of
        coefficients multiplying the basis functions defined by @a fluxelem. In
        other words, @a flux is the local vector for a FE space with
        dim*(dim+1)/2 vector components, based on the finite element @a fluxelem.
        The number of components, dim*(dim+1)/2 is such that it represents the
        symmetric part of the (symmetric) stress tensor. The order of the
        components is: s_xx, s_yy, s_xy in 2D, and s_xx, s_yy, s_zz, s_xy, s_xz,
        s_yz in 3D. */
    double ComputeFluxEnergy(const mfem::FiniteElement &fluxelem,
                             mfem::ElementTransformation &Trans,
                             mfem::Vector &flux,
                             mfem::Vector *d_energy = NULL) override;

    //! \brief destructor
    ~MGISIntegrator() override;

   protected:
    double q_lambda, q_mu;
    mfem::Coefficient *lambda, *mu;

   private:
#ifndef MFEM_THREAD_SAFE
    mfem::Vector shape;
    mfem::DenseMatrix dshape, gshape, pelmat;
    mfem::Vector divshape;
#endif

  }; // end of struct MGISIntegrator

}  // end of namesapce mfem_mgis

#endif /* LIB_MFEM_MGIS_MGISINTEGRATOR_HXX */
