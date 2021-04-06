/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BoundaryUtilities.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(const Parameters& p) {
    using SequentialImplementation =
        NonLinearEvolutionProblemImplementation<false>;
    const auto h =
        mgis::behaviour::fromString(get<std::string>(p, "Hypothesis"));
    const auto fed = std::make_shared<FiniteElementDiscretization>(
        extract(p, {"Parallel", "MeshFileName", "FiniteElementFamily",
                    "FiniteElementOrder", "UnknownsSize"}));
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      using ParallelImplementation =
          NonLinearEvolutionProblemImplementation<true>;
      this->pimpl = std::make_unique<ParallelImplementation>(fed, h, p);
#else
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      this->pimpl = std::make_unique<SequentialImplementation>(fed, h, p);
    }
  }  // end of NonLinearEvolutionProblem

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    using SequentialImplementation =
        NonLinearEvolutionProblemImplementation<false>;
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      using ParallelImplementation =
          NonLinearEvolutionProblemImplementation<true>;
      this->pimpl = std::make_unique<ParallelImplementation>(fed, h, p);
#else
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      this->pimpl = std::make_unique<SequentialImplementation>(fed, h, p);
    }
  }  // end of NonLinearEvolutionProblem

  FiniteElementDiscretization&
  NonLinearEvolutionProblem::getFiniteElementDiscretization() {
    return this->pimpl->getFiniteElementDiscretization();
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblem::getFiniteElementDiscretizationPointer() {
    return this->pimpl->getFiniteElementDiscretizationPointer();
  }  // end of getFiniteElementDiscretizationPointer

  void NonLinearEvolutionProblem::setSolverParameters(
      const Parameters& params) {
    return this->pimpl->setSolverParameters(params);
  }  // end of setSolverParameters

  void NonLinearEvolutionProblem::setLinearSolver(std::string_view n,
                                                  const Parameters& params) {
    this->pimpl->setLinearSolver(n, params);
  }  // end of setLinearSolver

  bool NonLinearEvolutionProblem::solve(const real t, const real dt) {
    this->setup(t, dt);
    return this->pimpl->solve(t, dt);
  }  // end of solve

  std::vector<size_type> NonLinearEvolutionProblem::getMaterialIdentifiers()
      const {
    return this->pimpl->getMaterialIdentifiers();
  }  // end of getMaterialIdentifiers

  void NonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->pimpl->addBoundaryCondition(std::move(bc));
  }  // end of NonLinearEvolutionProblem::addBoundaryCondition

  void NonLinearEvolutionProblem::addPostProcessing(
      const std::function<void(const real, const real)>& p) {
    this->pimpl->addPostProcessing(p);
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem::addPostProcessing(std::string_view n,
                                                    const Parameters& p) {
    this->pimpl->addPostProcessing(n, p);
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem::executePostProcessings(const real t,
                                                         const real dt) {
    this->pimpl->executePostProcessings(t, dt);
  }

  void NonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
                                                         const size_type l,
                                                         const std::string& m,
                                                         const std::string& b) {
    this->pimpl->addBehaviourIntegrator(n, l, m, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblem::getMaterial(
      const size_type m) const {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblem::getMaterial(const size_type m) {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  const BehaviourIntegrator& NonLinearEvolutionProblem::getBehaviourIntegrator(
      const size_type m) const {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  BehaviourIntegrator& NonLinearEvolutionProblem::getBehaviourIntegrator(
      const size_type m) {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblem::update() {
    this->pimpl->update();
  }  // end of update

  void NonLinearEvolutionProblem::revert() {
    this->pimpl->revert();
  }  // end of revert

  void NonLinearEvolutionProblem::setup(const real, const real) {
  }  // end of setup

  template <bool parallel>
  static NonLinearEvolutionProblemImplementation<parallel>&
  getImplementationInternal(AbstractNonLinearEvolutionProblem* const p) {
    auto* const pi =
        dynamic_cast<NonLinearEvolutionProblemImplementation<parallel>*>(p);
    if (pi == nullptr) {
      mgis::raise(
          "NonLinearEvolutionProblem::getImplementation: "
          "invalid call");
    }
    return *pi;
  }  // end of getImplementationInternal

  template <bool parallel>
  static const NonLinearEvolutionProblemImplementation<parallel>&
  getImplementationInternal(const AbstractNonLinearEvolutionProblem* const p) {
    const auto* const pi =
        dynamic_cast<const NonLinearEvolutionProblemImplementation<parallel>*>(
            p);
    if (pi == nullptr) {
      mgis::raise(
          "NonLinearEvolutionProblem::getImplementation: "
          "invalid call");
    }
    return *pi;
  }  // end of getImplementationInternal

  template <>
  const NonLinearEvolutionProblemImplementation<true>&
  NonLinearEvolutionProblem::getImplementation() const {
#ifdef MFEM_USE_MPI
    return getImplementationInternal<true>(this->pimpl.get());
#else  /* MFEM_USE_MPI */
    mgis::raise(
        "NonLinearEvolutionProblem::getImplementation: "
        "invalid call");
#endif /* MFEM_USE_MPI */
  }    // end of getImplementation

  template <>
  NonLinearEvolutionProblemImplementation<true>&
  NonLinearEvolutionProblem::getImplementation() {
#ifdef MFEM_USE_MPI
    return getImplementationInternal<true>(this->pimpl.get());
#else  /* MFEM_USE_MPI */
    mgis::raise(
        "NonLinearEvolutionProblem::getImplementation: "
        "invalid call");
#endif /* MFEM_USE_MPI */
  }    // end of getImplementation

  template <>
  const NonLinearEvolutionProblemImplementation<false>&
  NonLinearEvolutionProblem::getImplementation() const {
    return getImplementationInternal<false>(this->pimpl.get());
  }

  template <>
  NonLinearEvolutionProblemImplementation<false>&
  NonLinearEvolutionProblem::getImplementation() {
    return getImplementationInternal<false>(this->pimpl.get());
  }

  NonLinearEvolutionProblem::~NonLinearEvolutionProblem() = default;

  std::vector<std::pair<size_type, size_type>> buildFacesDescription(
      NonLinearEvolutionProblem& p, const size_type bid) {
    auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return buildFacesDescription(p.getImplementation<true>(), bid);
#else
      mgis::raise(
          "computeResultantForceOnBoundary: "
          "unsupported parallel computations");
#endif
    }
    return buildFacesDescription(p.getImplementation<false>(), bid);
  }  // end of buildFacesDescription

  void computeResultantForceOnBoundary(
      mfem::Vector& F,
      NonLinearEvolutionProblem& p,
      const std::vector<std::pair<size_type, size_type>>& faces) {
    auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      computeResultantForceOnBoundary(F, p.getImplementation<true>(), faces);
#else
      mgis::raise(
          "computeResultantForceOnBoundary: "
          "unsupported parallel computations");
#endif
    } else {
      computeResultantForceOnBoundary(F, p.getImplementation<false>(), faces);
    }
  }  // end of computeResultantForceOnBoundary

}  // end of namespace mfem_mgis