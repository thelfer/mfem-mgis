/*!
 * \file   src/UniformImposedPressureBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#include <iterator>
#include <iostream>

#include <algorithm>
#include "mfem/fem/nonlinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pnonlinearform.hpp"
#endif /* MFEM_USE_MPI */
#include "mfem/fem/nonlininteg.hpp"
#include "MFEMMGIS/Parameter.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"

namespace mfem_mgis {

  //! \brief nonlinear form implementing a uniform imposed pressure
  struct UniformImposedPressureBoundaryCondition::
      UniformImposedPressureNonlinearFormIntegratorBase
      : public NonlinearFormIntegrator {
    //! \brief constructor
    UniformImposedPressureNonlinearFormIntegratorBase()
        : pressure(0) {
    }  // end of UniformImposedPressureNonlinearFormIntegratorBase
    void setPressure(const real pr) { this->pressure = pr; }
    //! \brief destructor
    virtual ~UniformImposedPressureNonlinearFormIntegratorBase() = default;

   protected:
    real pressure;
  };

  //! \brief nonlinear form implementing a uniform imposed pressure
  struct UniformImposedPressureBoundaryCondition::
      UniformImposedPressureNonlinearFormIntegrator final
      : public UniformImposedPressureBoundaryCondition::
            UniformImposedPressureNonlinearFormIntegratorBase {
    //! \brief constructor
    UniformImposedPressureNonlinearFormIntegrator() {
    }  // end of UniformImposedPressureNonlinearFormIntegrator
    // MFEM API
    void AssembleElementVector(const mfem::FiniteElement &e,
                               mfem::ElementTransformation &tr,
                               const mfem::Vector &,
                               mfem::Vector &R) override {
      const int dim = e.GetDim() + 1;
      const int nnodes = e.GetDof();
      mfem::Vector n(dim);
      if (dim == 1) {
        n[0] = 1.0;
      }
#ifdef MFEM_THREAD_SAFE
      mfem::Vector shape;
#endif
      //
      shape.SetSize(nnodes);
      // initialize the residual
      R.SetSize(nnodes * dim);
      R = 0.0;
      // selection of the integration rule
      const auto *ir = this->IntRule;
      if (ir == nullptr) {
        const int o = e.GetOrder();
        ir = &mfem::IntRules.Get(e.GetGeomType(), o);
      }
      // loop over the integration point
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const auto &ip = ir->IntPoint(i);
        tr.SetIntPoint(&ip);
        // computation of the normal
        if (dim > 1) {
          CalcOrtho(tr.Jacobian(), n);
        }
        // computation of the external forces
        e.CalcShape(ip, shape);
        for (int ni = 0; ni != nnodes; ++ni) {
          for (int c = 0; c != dim; ++c) {
            R[ni + c * nnodes] +=
                ip.weight * (this->pressure) * n[c] * shape[ni];
          }
        }
      }
    }  // end of AssembleElementVector

    void AssembleElementGrad(const mfem::FiniteElement &e,
                             mfem::ElementTransformation &,
                             const mfem::Vector &,
                             mfem::DenseMatrix &K) override {
      const int dim = e.GetDim() + 1;
      const int nnodes = e.GetDof();
      K.SetSize(nnodes * dim, nnodes * dim);
      K = 0.0;
    }
    //! \brief destructor
    virtual ~UniformImposedPressureNonlinearFormIntegrator() = default;

#ifndef MFEM_THREAD_SAFE
   private:
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
#endif
  };  // end of UniformImposedPressureNonlinearFormIntegrator

  UniformImposedPressureBoundaryCondition::
      UniformImposedPressureBoundaryCondition(
          AbstractNonLinearEvolutionProblem &p, const Parameters &params)
      : finiteElementDiscretization(p.getFiniteElementDiscretizationPointer()),
        bids(getBoundariesIdentifiers(p, params, false)),
        prfct(get<std::function<real(const real)>>(params, "LoadingEvolution")),
        nfi(new UniformImposedPressureNonlinearFormIntegrator) {}

  UniformImposedPressureBoundaryCondition::
      UniformImposedPressureBoundaryCondition(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const size_type bid,
          std::function<real(const real)> prvalues)
      : finiteElementDiscretization(fed),
        bids(1, bid),
        prfct(prvalues),
        nfi(new UniformImposedPressureNonlinearFormIntegrator) {}

  UniformImposedPressureBoundaryCondition::
      UniformImposedPressureBoundaryCondition(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const std::string_view bid,
          std::function<real(const real)> prvalues)
      : finiteElementDiscretization(fed),
        bids(fed->getBoundariesIdentifiers(bid)),
        prfct(prvalues),
        nfi(new UniformImposedPressureNonlinearFormIntegrator) {}

  void UniformImposedPressureBoundaryCondition::setup(const real t,
                                                      const real dt) {
    this->nfi->setPressure(this->prfct(t + dt));
  }

#ifdef MFEM_USE_MPI
  void UniformImposedPressureBoundaryCondition::addNonlinearFormIntegrator(
      NonlinearForm<true> &f) {
    auto &m = this->finiteElementDiscretization->getMesh<true>();
    this->boundaries_markers =
        mfem::Array<mfem_mgis::size_type>(m.bdr_attributes.Max());
    this->boundaries_markers = 0;
    for (const auto &bid : bids) {
      this->boundaries_markers[bid - 1] = 1;
    }
    f.AddBoundaryIntegrator(this->nfi, this->boundaries_markers);
    this->shallFreeIntegrator = false;
  }    // end of addNonlinearFormIntegrator
#endif /* MFEM_USE_MPI */

  void UniformImposedPressureBoundaryCondition::addNonlinearFormIntegrator(
      NonlinearForm<false> &f) {
    auto &m = this->finiteElementDiscretization->getMesh<false>();
    this->boundaries_markers =
        mfem::Array<mfem_mgis::size_type>(m.bdr_attributes.Max());
    this->boundaries_markers = 0;
    for (const auto &bid : bids) {
      this->boundaries_markers[bid - 1] = 1;
    }
    f.AddBoundaryIntegrator(this->nfi, this->boundaries_markers);
    this->shallFreeIntegrator = false;
  }  // end of addNonlinearFormIntegrator

  UniformImposedPressureBoundaryCondition::
      ~UniformImposedPressureBoundaryCondition() {
    if (this->shallFreeIntegrator) {
      std::free(this->nfi);
    }
  }

}  // end of namespace mfem_mgis