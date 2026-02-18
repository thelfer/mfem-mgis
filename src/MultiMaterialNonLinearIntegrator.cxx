/*!
 * \file   src/MultiMaterialNonLinearIntegrator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "mfem/fem/linearform.hpp"
#include "mfem/fem/bilinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/mesh/pmesh.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple test
   * \param[in] i: pointer to behaviour integrator
   * \param[in] n: name of the calling method
   * \param[in] m: material id
   */
  static void checkIfBehaviourIntegratorIsDefined(
      const AbstractBehaviourIntegrator* const i,
      const char* const n,
      const size_type m) {
    if (i == nullptr) {
      raise("MultiMaterialNonLinearIntegrator::" + std::string(n) +
            ": no behaviour integrator associated with material '" +
            std::to_string(m) + "'");
    }
  }  // end if checkIfBehaviourIntegratorIsDefined

  struct InternalForcesIntegrator final : public LinearFormIntegrator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    InternalForcesIntegrator(
        const FiniteElementSpace<true>& fes,
        const mfem::Vector& u,
        std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>& bis)
        : pfespace(&fes), unknowns(u), behaviour_integrators(bis) {}
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    InternalForcesIntegrator(
        const FiniteElementSpace<false>& fes,
        const mfem::Vector& u,
        std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>& bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& F) override {
      const auto m = tr.Attribute;
      const auto& bi = this->behaviour_integrators[m];
      checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementVector", m);
      if (bi->requiresCurrentSolutionForResidualAssembly()) {
        auto vdofs = mfem::Array<int>{};
#ifdef MFEM_USE_MPI
        if (pfespace != nullptr) {
          pfespace->GetElementVDofs(tr.ElementNo, vdofs);
        } else {
          fespace->GetElementVDofs(tr.ElementNo, vdofs);
        }
#else  /* MFEM_USE_MPI */
        fespace->GetElementVDofs(tr.ElementNo, vdofs);
#endif /* MFEM_USE_MPI */
        this->unknowns.GetSubVector(vdofs, this->Ue);
      }
      bi->updateResidual(F, e, tr, this->Ue);
    }  // end of AssembleRHSElementVect

   private:
#ifdef MFEM_USE_MPI
    //! \brief pointer to the finite element space for parallel resolutions
    const FiniteElementSpace<true>* const pfespace = nullptr;
#endif /* MFEM_USE_MPI */
    //! \brief pointer to the finite element space for sequential resolutions
    const FiniteElementSpace<false>* const fespace = nullptr;
    //! \brief unknowns
    const mfem::Vector& unknowns;
    //! \brief temporary vector containing the unknown of the current element
    mfem::Vector Ue;
    //! \brief list of behaviour integrators
    std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>&
        behaviour_integrators;
  };

  struct StiffnessMatrixIntegrator final : public BilinearFormIntegrator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    StiffnessMatrixIntegrator(
        const FiniteElementSpace<true>& fes,
        const mfem::Vector& u,
        std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>& bis)
        : pfespace(&fes), unknowns(u), behaviour_integrators(bis) {}
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] u: unknown
     * \param[in] fes: finite element space
     * \param[in] bis: behaviour integrators
     */
    StiffnessMatrixIntegrator(
        const FiniteElementSpace<false>& fes,
        const mfem::Vector& u,
        std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>& bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleElementMatrix(const mfem::FiniteElement& e,
                               mfem::ElementTransformation& tr,
                               mfem::DenseMatrix& K) override {
      const auto m = tr.Attribute;
      const auto& bi = this->behaviour_integrators[m];
      checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementGrad", m);
      if (bi->requiresCurrentSolutionForJacobianAssembly()) {
        auto vdofs = mfem::Array<int>{};
#ifdef MFEM_USE_MPI
        if (pfespace != nullptr) {
          pfespace->GetElementVDofs(tr.ElementNo, vdofs);
        } else {
          fespace->GetElementVDofs(tr.ElementNo, vdofs);
        }
#else  /* MFEM_USE_MPI */
        fespace->GetElementVDofs(tr.ElementNo, vdofs);
#endif /* MFEM_USE_MPI */
        this->unknowns.GetSubVector(vdofs, this->Ue);
      }
      bi->updateJacobian(K, e, tr, this->Ue);
    }  // end of AssembleElementMatrix

   private:
#ifdef MFEM_USE_MPI
    //! \brief pointer to the finite element space for parallel resolutions
    const FiniteElementSpace<true>* const pfespace = nullptr;
#endif /* MFEM_USE_MPI */
    //! \brief pointer to the finite element space for sequential resolutions
    const FiniteElementSpace<false>* const fespace = nullptr;
    //! \brief unknowns
    const mfem::Vector& unknowns;
    //! \brief temporary vector containing the unknown of the current element
    mfem::Vector Ue;
    //! \brief list of behaviour integrators
    std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>&
        behaviour_integrators;
  };  // end of StiffnessMatrixIntegrator

  MultiMaterialNonLinearIntegrator::MultiMaterialNonLinearIntegrator(
      std::shared_ptr<const FiniteElementDiscretization> fed,
      const Hypothesis h)
      : fe_discretization(fed), hypothesis(h) {
    if (this->fe_discretization->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = this->fe_discretization->getMesh<true>();
      // shifting by one allows to directly use the material id to get
      // the behaviour integrators in all the other methods.
      // However, behaviour_integrators[0] will always be null.
      this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto& mesh = this->fe_discretization->getMesh<false>();
      // shifting by one allows to directly use the material id to get
      // the behaviour integrators in all the other methods.
      // However, behaviour_integrators[0] will always be null.
      this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
    }
  }  // end of MultiMaterialNonLinearIntegrator

  bool MultiMaterialNonLinearIntegrator::integrate(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      const IntegrationType it) {
    const auto m = tr.Attribute;
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "integrate", m);
    return bi->integrate(e, tr, U, it);
  }  // end of integrate

  void MultiMaterialNonLinearIntegrator::AssembleElementVector(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::Vector& F) {
    const auto m = tr.Attribute;
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementVector", m);
    if (usePETSc()) {
      MFEM_VERIFY(bi->integrate(
                      e, tr, U,
                      IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR),
                  "ERROR Behaviour");
    }
    bi->updateResidual(F, e, tr, U);
  }  // end of AssembleElementVector

  void MultiMaterialNonLinearIntegrator::AssembleElementGrad(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::DenseMatrix& K) {
    const auto m = tr.Attribute;
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "AssembleElementGrad", m);
    bi->updateJacobian(K, e, tr, U);
  }  // end of AssembleElementGrad

  void MultiMaterialNonLinearIntegrator::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    if (this->behaviour_integrators[m] != nullptr) {
      raise(
          "MultiMaterialNonLinearIntegrator::addBehaviourIntegrator: "
          "integrator already defined for material '" +
          std::to_string(m) + "'");
    }
    const auto& f = BehaviourIntegratorFactory::get(this->hypothesis);
    this->behaviour_integrators[m] =
        f.generate(n, *(this->fe_discretization), m,
                   mfem_mgis::load(l, b, this->hypothesis));
  }  // end of addBehaviourIntegrator

  const Material& MultiMaterialNonLinearIntegrator::getMaterial(
      const size_type m) const {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
    return bi->getMaterial();
  }  // end of getMaterial

  Material& MultiMaterialNonLinearIntegrator::getMaterial(const size_type m) {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getMaterial", m);
    return bi->getMaterial();
  }  // end of getMaterial

  const AbstractBehaviourIntegrator&
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      const size_type m) const {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getBehaviourIntegrator", m);
    return *bi;
  }  // end of getBehaviourIntegrator

  AbstractBehaviourIntegrator&
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(const size_type m) {
    const auto& bi = this->behaviour_integrators[m];
    checkIfBehaviourIntegratorIsDefined(bi.get(), "getBehaviourIntegrator", m);
    return *bi;
  }  // end of getBehaviourIntegrator

  void MultiMaterialNonLinearIntegrator::setTimeIncrement(const real dt) {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->setTimeIncrement(dt);
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setMacroscopicGradients(
      std::span<const real> g) {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->setMacroscopicGradients(g);
      }
    }
  }  // end of setMacroscopicGradients

  void MultiMaterialNonLinearIntegrator::setup(const real t, const real dt) {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->setup(t, dt);
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::revert() {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->revert();
      }
    }
  }  // end of revert

  void MultiMaterialNonLinearIntegrator::update() {
    for (auto& bi : this->behaviour_integrators) {
      if (bi != nullptr) {
        bi->update();
      }
    }
  }  // end of update

  std::vector<size_type>
  MultiMaterialNonLinearIntegrator::getAssignedMaterialsIdentifiers() const {
    std::vector<size_type> mids;
    for (size_type i = 0; i != this->behaviour_integrators.size(); ++i) {
      if (this->behaviour_integrators[i] != nullptr) {
        mids.push_back(i);
      }
    }
    return mids;
  }  // end of getAssignedMaterialsIdentifiers

  LinearizedOperators MultiMaterialNonLinearIntegrator::getLinearizedOperators(
      const mfem::Vector& u) {
    if (this->fe_discretization->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& fespace =
          this->fe_discretization->getFiniteElementSpace<true>();
      return LinearizedOperators{
          .K = std::make_unique<StiffnessMatrixIntegrator>(
              fespace, u, this->behaviour_integrators),
          .mFi = std::make_unique<InternalForcesIntegrator>(
              fespace, u, this->behaviour_integrators)};
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& fespace =
        this->fe_discretization->getFiniteElementSpace<false>();
    return LinearizedOperators{
        .K = std::make_unique<StiffnessMatrixIntegrator>(
            fespace, u, this->behaviour_integrators),
        .mFi = std::make_unique<InternalForcesIntegrator>(
            fespace, u, this->behaviour_integrators)};
  }  // end of getLinearizedOperators

  MultiMaterialNonLinearIntegrator::~MultiMaterialNonLinearIntegrator() =
      default;

}  // end of namespace mfem_mgis
