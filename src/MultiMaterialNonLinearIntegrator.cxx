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
  static void checkIfBehaviourIntegratorsAreDefined(
      const std::vector<
          std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>& bis,
      const char* const n,
      const size_type m) {
    if (m > bis.size()) {
      raise("MultiMaterialNonLinearIntegrator::" + std::string(n) +
            ": no behaviour integrator associated with material '" +
            std::to_string(m) + "'");
    }
    // ok this is paranoïac, but does not hurt
    if (bis.at(m).front().get() == nullptr) {
      raise("MultiMaterialNonLinearIntegrator::" + std::string(n) +
            ": invalid behaviour integrator associated with material '" +
            std::to_string(m) + "'");
    }
  }  // end if checkIfBehaviourIntegratorsAreDefined

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
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
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
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleRHSElementVect(const mfem::FiniteElement& e,
                                mfem::ElementTransformation& tr,
                                mfem::Vector& F) override {
      const auto m = tr.Attribute;
      checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                            "AssembleElementVector", m);
      for (auto& bi : this->behaviour_integrators.at(m)) {
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
      }
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
    std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
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
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
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
        std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
            bis)
        : fespace(&fes), unknowns(u), behaviour_integrators(bis) {}
    //
    void AssembleElementMatrix(const mfem::FiniteElement& e,
                               mfem::ElementTransformation& tr,
                               mfem::DenseMatrix& K) override {
      const auto m = tr.Attribute;
      checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                            "AssembleElementGrad", m);
      for (auto& bi : this->behaviour_integrators.at(m)) {
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
      }
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
    std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>&
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
      // However, behaviour_integrators[0] will always be empty.
      this->behaviour_integrators.resize(mesh.attributes.Max() + 1);
    }
  }  // end of MultiMaterialNonLinearIntegrator

  bool MultiMaterialNonLinearIntegrator::integrate(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      const IntegrationType it) {
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "integrate", m);
    for (const auto& bi : this->behaviour_integrators.at(m)) {
      if (!bi->integrate(e, tr, U, it)) {
        return false;
      }
    }
    return true;
  }  // end of integrate

  void MultiMaterialNonLinearIntegrator::AssembleElementVector(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::Vector& F) {
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "AssembleElementVector", m);
    for (const auto& bi : this->behaviour_integrators.at(m)) {
      if (usePETSc()) {
        MFEM_VERIFY(
            bi->integrate(
                e, tr, U,
                IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR),
            "ERROR Behaviour");
      }
      bi->updateResidual(F, e, tr, U);
    }
  }  // end of AssembleElementVector

  void MultiMaterialNonLinearIntegrator::AssembleElementGrad(
      const mfem::FiniteElement& e,
      mfem::ElementTransformation& tr,
      const mfem::Vector& U,
      mfem::DenseMatrix& K) {
    const auto m = tr.Attribute;
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "AssembleElementGrad", m);
    for (const auto& bi : this->behaviour_integrators.at(m)) {
      bi->updateJacobian(K, e, tr, U);
    }
  }  // end of AssembleElementGrad

  size_type MultiMaterialNonLinearIntegrator::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    const auto& f = BehaviourIntegratorFactory::get(this->hypothesis);
    auto& bis = this->behaviour_integrators[m];
    const auto s = static_cast<size_type>(bis.size());
    bis.push_back(f.generate(n, *(this->fe_discretization), m,
                             mfem_mgis::load(l, b, this->hypothesis)));
    return s;
  }  // end of addBehaviourIntegrator

  OptionalReference<const Material>
  MultiMaterialNonLinearIntegrator::getMaterial(
      Context& ctx, const size_type m, const size_type b) const noexcept {
    const auto obi = this->getBehaviourIntegrator(ctx, m, b);
    if (isInvalid(obi)) {
      return {};
    }
    return {&(obi->getMaterial())};
  }  // end of getMaterial

  OptionalReference<Material> MultiMaterialNonLinearIntegrator::getMaterial(
      Context& ctx, const size_type m, const size_type b) noexcept {
    auto obi = this->getBehaviourIntegrator(ctx, m, b);
    if (isInvalid(obi)) {
      return {};
    }
    return {&(obi->getMaterial())};
  }  // end of getMaterial

  OptionalReference<const AbstractBehaviourIntegrator>
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      Context& ctx, const size_type m, const size_type b) const noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    const auto& bis = this->behaviour_integrators.at(m);
    if (b > bis.size()) {
      return ctx.registerErrorMessage("invalid behaviour index '" +
                                      std::to_string(b) + "'");
    }
    return {bis.at(b).get()};
  }  // end of getBehaviourIntegrator

  OptionalReference<AbstractBehaviourIntegrator>
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      Context& ctx, const size_type m, const size_type b) noexcept {
    if ((m == 0) || (m >= this->behaviour_integrators.size())) {
      return ctx.registerErrorMessage("invalid material index '" +
                                      std::to_string(m) + "'");
    }
    auto& bis = this->behaviour_integrators.at(m);
    if (b > bis.size()) {
      return ctx.registerErrorMessage("invalid behaviour index '" +
                                      std::to_string(b) + "'");
    }
    return {bis.at(b).get()};
  }  // end of getBehaviourIntegrator

  const Material& MultiMaterialNonLinearIntegrator::getMaterial(
      const size_type m) const {
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "getMaterial", m);
    const auto& bis = this->behaviour_integrators[m];
    return bis.front()->getMaterial();
  }  // end of getMaterial

  Material& MultiMaterialNonLinearIntegrator::getMaterial(const size_type m) {
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "getMaterial", m);
    const auto& bis = this->behaviour_integrators[m];
    return bis.front()->getMaterial();
  }  // end of getMaterial

  const AbstractBehaviourIntegrator&
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(
      const size_type m) const {
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "getBehaviourIntegrator", m);
    const auto& bis = this->behaviour_integrators[m];
    return *(bis.front());
  }  // end of getBehaviourIntegrator

  AbstractBehaviourIntegrator&
  MultiMaterialNonLinearIntegrator::getBehaviourIntegrator(const size_type m) {
    checkIfBehaviourIntegratorsAreDefined(this->behaviour_integrators,
                                          "getBehaviourIntegrator", m);
    const auto& bis = this->behaviour_integrators[m];
    return *(bis.front());
  }  // end of getBehaviourIntegrator

  real MultiMaterialNonLinearIntegrator::getTimeIncrement() const noexcept {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
	if (bi != nullptr) {
	  return bi->getTimeIncrement();
	}
      }
    }
    return real{};
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setTimeIncrement(const real dt) {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        if (bi != nullptr) {
          bi->setTimeIncrement(dt);
        }
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::setMacroscopicGradients(
      std::span<const real> g) {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->setMacroscopicGradients(g);
      }
    }
  }  // end of setMacroscopicGradients

  void MultiMaterialNonLinearIntegrator::setup(const real t, const real dt) {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->setup(t, dt);
      }
    }
  }  // end of setTimeIncrement

  void MultiMaterialNonLinearIntegrator::revert() {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->revert();
      }
    }
  }  // end of revert

  void MultiMaterialNonLinearIntegrator::update() {
    for (auto& bis : this->behaviour_integrators) {
      for (auto& bi : bis) {
        bi->update();
      }
    }
  }  // end of update

  std::vector<size_type>
  MultiMaterialNonLinearIntegrator::getAssignedMaterialsIdentifiers() const {
    std::vector<size_type> mids;
    for (size_type i = 0; i != this->behaviour_integrators.size(); ++i) {
      if (!this->behaviour_integrators[i].empty()) {
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
