/*!
 * \file   src/UniformHeatSourceBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#include <iterator>
#include <iostream>

#include <algorithm>
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/nonlinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/plinearform.hpp"
#include "mfem/fem/pnonlinearform.hpp"
#endif /* MFEM_USE_MPI */
#include "mfem/fem/nonlininteg.hpp"
#include "MFEMMGIS/Parameter.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformHeatSourceBoundaryCondition.hxx"

namespace mfem_mgis {

  //! \brief nonlinear form implementing a uniform imposed heat source
  struct UniformHeatSourceBoundaryCondition::
      UniformHeatSourceFormIntegratorBase {
    //! \brief constructor
    UniformHeatSourceFormIntegratorBase()
        : heat_source(0) {
    }  // end of UniformHeatSourceFormIntegratorBase
    void setHeatSource(const real pr) { this->heat_source = pr; }
    //! \brief destructor
    virtual ~UniformHeatSourceFormIntegratorBase() = default;

   protected:
    //
    [[nodiscard]] virtual const mfem::IntegrationRule *getIntegrationRule(
        const mfem::FiniteElement &e) const noexcept = 0;
    //
    void computeResidual(mfem::Vector &R,
                         const mfem::FiniteElement &e,
                         mfem::ElementTransformation &tr) {
      const int nnodes = e.GetDof();
#ifdef MFEM_THREAD_SAFE
      mfem::Vector shape;
#endif
      //
      shape.SetSize(nnodes);
      // initialize the residual
      R.SetSize(nnodes);
      R = 0.0;
      // selection of the integration rule
      const auto *ir = this->getIntegrationRule(e);
      // loop over the integration point
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const auto &ip = ir->IntPoint(i);
        tr.SetIntPoint(&ip);
        // computation of the external forces
        e.CalcShape(ip, shape);
        for (int ni = 0; ni != nnodes; ++ni) {
          R[ni] = ip.weight * tr.Weight() * (this->heat_source) * shape[ni];
        }
      }
    }  // end of computeResidual

    real heat_source;

#ifndef MFEM_THREAD_SAFE
   private:
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
#endif
  };

  //! \brief linear form implementing a uniform imposed heat source
  struct UniformHeatSourceBoundaryCondition::
      UniformHeatSourceLinearFormIntegrator final
      : public LinearFormIntegrator,
        public UniformHeatSourceBoundaryCondition::
            UniformHeatSourceFormIntegratorBase {
    //! \brief constructor
    UniformHeatSourceLinearFormIntegrator() {
    }  // end of UniformHeatSourceLinearFormIntegrator
    // MFEM API
    void AssembleRHSElementVect(const mfem::FiniteElement &e,
                                mfem::ElementTransformation &tr,
                                mfem::Vector &R) override {
      this->computeResidual(R, e, tr);
    }  // end of AssembleElementVector
    //! \brief destructor
    virtual ~UniformHeatSourceLinearFormIntegrator() = default;

   protected:
    [[nodiscard]] const mfem::IntegrationRule *getIntegrationRule(
        const mfem::FiniteElement &e) const noexcept override {
      const auto *ir = this->IntRule;
      if (ir == nullptr) {
        const int o = e.GetOrder();
        return &mfem::IntRules.Get(e.GetGeomType(), o);
      }
      return ir;
    }  // end of getIntegrationRule
  };   // end of UniformHeatSourceLinearFormIntegrator

  //! \brief nonlinear form implementing a uniform imposed heat source
  struct UniformHeatSourceBoundaryCondition::
      UniformHeatSourceNonlinearFormIntegrator final
      : public NonlinearFormIntegrator,
        public UniformHeatSourceBoundaryCondition::
            UniformHeatSourceFormIntegratorBase {
    //! \brief constructor
    UniformHeatSourceNonlinearFormIntegrator() {
    }  // end of UniformHeatSourceNonlinearFormIntegrator
    // MFEM API
    void AssembleElementVector(const mfem::FiniteElement &e,
                               mfem::ElementTransformation &tr,
                               const mfem::Vector &,
                               mfem::Vector &R) override {
      this->computeResidual(R, e, tr);
    }  // end of AssembleElementVector

    void AssembleElementGrad(const mfem::FiniteElement &e,
                             mfem::ElementTransformation &,
                             const mfem::Vector &,
                             mfem::DenseMatrix &K) override {
      const int nnodes = e.GetDof();
      K.SetSize(nnodes, nnodes);
      K = 0.0;
    }
    //! \brief destructor
    virtual ~UniformHeatSourceNonlinearFormIntegrator() = default;

   protected:
    [[nodiscard]] const mfem::IntegrationRule *getIntegrationRule(
        const mfem::FiniteElement &e) const noexcept override {
      const auto *ir = this->IntRule;
      if (ir == nullptr) {
        const int o = e.GetOrder();
        return &mfem::IntRules.Get(e.GetGeomType(), o);
      }
      return ir;
    }  // end of getIntegrationRule
  };  // end of UniformHeatSourceNonlinearFormIntegrator

  UniformHeatSourceBoundaryCondition::UniformHeatSourceBoundaryCondition(
      AbstractNonLinearEvolutionProblem &p, const Parameters &params)
      : finiteElementDiscretization(p.getFiniteElementDiscretizationPointer()),
        mids(getBoundariesIdentifiers(p, params, false)),
        qfct(get<std::function<real(const real)>>(params, "LoadingEvolution")),
        nfi(new UniformHeatSourceNonlinearFormIntegrator) {}

  UniformHeatSourceBoundaryCondition::UniformHeatSourceBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type mid,
      std::function<real(const real)> prvalues)
      : finiteElementDiscretization(fed),
        mids(1, mid),
        qfct(prvalues),
        nfi(new UniformHeatSourceNonlinearFormIntegrator) {}

  UniformHeatSourceBoundaryCondition::UniformHeatSourceBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const std::string_view mid,
      std::function<real(const real)> prvalues)
      : finiteElementDiscretization(fed),
        mids(fed->getMaterialsIdentifiers(mid)),
        qfct(prvalues),
        nfi(new UniformHeatSourceNonlinearFormIntegrator) {}

  void UniformHeatSourceBoundaryCondition::setup(const real t, const real dt) {
    this->nfi->setHeatSource(this->qfct(t + dt));
  }

#ifdef MFEM_USE_MPI
  void UniformHeatSourceBoundaryCondition::addNonlinearFormIntegrator(
      NonlinearForm<true> &f) {
    auto &m = this->finiteElementDiscretization->getMesh<true>();
    this->materials_markers =
        mfem::Array<mfem_mgis::size_type>(m.attributes.Max());
    this->materials_markers = 0;
    for (const auto &mid : mids) {
      this->materials_markers[mid - 1] = 1;
    }
    f.AddDomainIntegrator(this->nfi, this->materials_markers);
    this->shallFreeIntegrator = false;
  }    // end of addNonlinearFormIntegrator
#endif /* MFEM_USE_MPI */

  void UniformHeatSourceBoundaryCondition::addNonlinearFormIntegrator(
      NonlinearForm<false> &f) {
    auto &m = this->finiteElementDiscretization->getMesh<false>();
    this->materials_markers =
        mfem::Array<mfem_mgis::size_type>(m.attributes.Max());
    this->materials_markers = 0;
    for (const auto &mid : mids) {
      this->materials_markers[mid - 1] = 1;
    }
    f.AddDomainIntegrator(this->nfi, this->materials_markers);
    this->shallFreeIntegrator = false;
  }  // end of addNonlinearFormIntegrator

#ifdef MFEM_USE_MPI
  void UniformHeatSourceBoundaryCondition::addLinearFormIntegrator(
      LinearForm<true> &f, const real t, const real dt) {
    auto &m = this->finiteElementDiscretization->getMesh<true>();
    this->materials_markers =
        mfem::Array<mfem_mgis::size_type>(m.attributes.Max());
    this->materials_markers = 0;
    for (const auto &bid : mids) {
      this->materials_markers[bid - 1] = 1;
    }
    auto *form = new UniformHeatSourceLinearFormIntegrator();
    form->setHeatSource(this->qfct(t + dt));
    f.AddBoundaryIntegrator(form, this->materials_markers);
  }    // end of addLinearFormIntegrator
#endif /* MFEM_USE_MPI */

  void UniformHeatSourceBoundaryCondition::addLinearFormIntegrator(
      LinearForm<false> &f, const real t, const real dt) {
    auto &m = this->finiteElementDiscretization->getMesh<false>();
    this->materials_markers =
        mfem::Array<mfem_mgis::size_type>(m.attributes.Max());
    this->materials_markers = 0;
    for (const auto &bid : mids) {
      this->materials_markers[bid - 1] = 1;
    }
    auto *form = new UniformHeatSourceLinearFormIntegrator();
    form->setHeatSource(this->qfct(t + dt));
    f.AddBoundaryIntegrator(form, this->materials_markers);
  }  // end of addLinearFormIntegrator

  UniformHeatSourceBoundaryCondition::~UniformHeatSourceBoundaryCondition() {
    if (this->shallFreeIntegrator) {
      std::free(this->nfi);
    }
  }

}  // end of namespace mfem_mgis
