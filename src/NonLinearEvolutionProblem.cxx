/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/BoundaryUtilities.hxx"
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/DirichletBoundaryCondition.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  const char* const NonLinearEvolutionProblem::HypothesisParameter =
      "Hypothesis";

  std::vector<std::string> NonLinearEvolutionProblem::getParametersList() {
    auto params = FiniteElementDiscretization::getParametersList();
    auto params2 =
        NonLinearEvolutionProblemImplementationBase::getParametersList();
    params.insert(params.end(), params2.begin(), params2.end());
    params.push_back(NonLinearEvolutionProblem::HypothesisParameter);
    return params;
  }  // end of getParametersList

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(const Parameters& p) {
    using SequentialImplementation =
        NonLinearEvolutionProblemImplementation<false>;
    const auto h = mgis::behaviour::fromString(
        get<std::string>(p, NonLinearEvolutionProblem::HypothesisParameter));
    const auto fed = std::make_shared<FiniteElementDiscretization>(
        extract(p, FiniteElementDiscretization::getParametersList()));
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      using ParallelImplementation =
          NonLinearEvolutionProblemImplementation<true>;
      this->pimpl = std::make_unique<ParallelImplementation>(fed, h, p);
#else
      raise(
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
      raise(
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

  const FiniteElementDiscretization&
  NonLinearEvolutionProblem::getFiniteElementDiscretization() const {
    return this->pimpl->getFiniteElementDiscretization();
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblem::getFiniteElementDiscretizationPointer() {
    return this->pimpl->getFiniteElementDiscretizationPointer();
  }  // end of getFiniteElementDiscretizationPointer

  mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() {
    return this->pimpl->getUnknownsAtBeginningOfTheTimeStep();
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblem::getUnknownsAtBeginningOfTheTimeStep() const {
    return this->pimpl->getUnknownsAtBeginningOfTheTimeStep();
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep() {
    return this->pimpl->getUnknownsAtEndOfTheTimeStep();
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblem::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->pimpl->getUnknownsAtEndOfTheTimeStep();
  }  // end of getUnknownsAtEndOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblem::getUnknowns(
      const TimeStepStage ts) noexcept {
    return this->pimpl->getUnknowns(ts);
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblem::getUnknowns(
      const TimeStepStage ts) const noexcept {
    return this->pimpl->getUnknowns(ts);
  }  // end of getUnknownsAtEndOfTheTimeStep

  void NonLinearEvolutionProblem::setSolverParameters(
      const Parameters& params) {
    return this->pimpl->setSolverParameters(params);
  }  // end of setSolverParameters

  bool NonLinearEvolutionProblem::setLinearSolver(
      Context& ctx, LinearSolverHandler s) noexcept {
    return this->pimpl->setLinearSolver(ctx, std::move(s));
  }  // end of setLinearSolver

  bool
  NonLinearEvolutionProblem::areStiffnessOperatorsFromLastIterationAvailable()
      const noexcept {
    return this->pimpl->areStiffnessOperatorsFromLastIterationAvailable();
  }  // end of areStiffnessOperatorsFromLastIterationAvailable

  void NonLinearEvolutionProblem::setLinearSolver(std::string_view n,
                                                  const Parameters& params) {
    this->pimpl->setLinearSolver(n, params);
  }  // end of setLinearSolver

  void NonLinearEvolutionProblem::setPredictionPolicy(
      const PredictionPolicy& p) noexcept {
    this->pimpl->setPredictionPolicy(p);
  }  // end of setPredictionPolicy

  PredictionPolicy NonLinearEvolutionProblem::getPredictionPolicy()
      const noexcept {
    return this->pimpl->getPredictionPolicy();
  }  // end of getPredictionPolicy

  const std::vector<std::unique_ptr<DirichletBoundaryCondition>>&
  NonLinearEvolutionProblem::getDirichletBoundaryConditions() const noexcept {
    return this->pimpl->getDirichletBoundaryConditions();
  }  // end of getDirichletBoundaryConditions

  const std::vector<std::unique_ptr<AbstractBoundaryCondition>>&
  NonLinearEvolutionProblem::getBoundaryConditions() const noexcept {
    return this->pimpl->getBoundaryConditions();
  }  // end of getBoundaryConditions

  NonLinearResolutionOutput NonLinearEvolutionProblem::solve(const real t,
                                                             const real dt) {
    CatchTimeSection("NLEP::solve");
    this->setup(t, dt);
    return this->pimpl->solve(t, dt);
  }  // end of solve

  bool NonLinearEvolutionProblem::integrate(const mfem::Vector& U,
                                            const IntegrationType it,
                                            const std::optional<real> odt) {
    CatchTimeSection("NLEP::integrate");
    return this->pimpl->integrate(U, it, odt);
  }  // end of solve

  std::vector<size_type>
  NonLinearEvolutionProblem::getEssentialDegreesOfFreedom() const {
    return this->pimpl->getEssentialDegreesOfFreedom();
  }  // end of getEssentialDegreesOfFreedom

  std::optional<LinearizedOperators>
  NonLinearEvolutionProblem::getLinearizedOperators(
      Context& ctx, const mfem::Vector& U) noexcept {
    return this->pimpl->getLinearizedOperators(ctx, U);
  }  // end of getLinearizedOperators

  bool NonLinearEvolutionProblem::setMaterialsNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->pimpl->setMaterialsNames(ctx, ids);
  }  // end of setMaterialsNames

  bool NonLinearEvolutionProblem::setBoundariesNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->pimpl->setBoundariesNames(ctx, ids);
  }  // end of setBoundariesNames

  void NonLinearEvolutionProblem::setMaterialsNames(
      const std::map<size_type, std::string>& ids) {
    this->pimpl->setMaterialsNames(ids);
  }

  void NonLinearEvolutionProblem::setBoundariesNames(
      const std::map<size_type, std::string>& ids) {
    this->pimpl->setBoundariesNames(ids);
  }

  std::optional<size_type> NonLinearEvolutionProblem::getMaterialIdentifier(
      Context& ctx, const Parameter& m) const noexcept {
    return this->pimpl->getMaterialIdentifier(ctx, m);
  }

  std::optional<size_type> NonLinearEvolutionProblem::getBoundaryIdentifier(
      Context& ctx, const Parameter& b) const noexcept {
    return this->pimpl->getBoundaryIdentifier(ctx, b);
  }

  std::optional<std::vector<size_type>>
  NonLinearEvolutionProblem::getMaterialsIdentifiers(
      Context& ctx, const Parameter& m) const noexcept {
    return this->pimpl->getMaterialsIdentifiers(ctx, m);
  }

  std::optional<std::vector<size_type>>
  NonLinearEvolutionProblem::getBoundariesIdentifiers(
      Context& ctx, const Parameter& b) const noexcept {
    return this->pimpl->getBoundariesIdentifiers(ctx, b);
  }

  size_type NonLinearEvolutionProblem::getMaterialIdentifier(
      const Parameter& p) const {
    return this->pimpl->getMaterialIdentifier(p);
  }  // end of getMaterialIdentifier

  size_type NonLinearEvolutionProblem::getBoundaryIdentifier(
      const Parameter& p) const {
    return this->pimpl->getBoundaryIdentifier(p);
  }  // end of getBoundariesIdentifier

  std::vector<size_type> NonLinearEvolutionProblem::getMaterialsIdentifiers(
      const Parameter& p) const {
    return this->pimpl->getMaterialsIdentifiers(p);
  }  // end of getMaterialsIdentifiers

  std::vector<size_type> NonLinearEvolutionProblem::getBoundariesIdentifiers(
      const Parameter& p) const {
    return this->pimpl->getBoundariesIdentifiers(p);
  }  // end of getBoundariesIdentifiers

  std::vector<size_type>
  NonLinearEvolutionProblem::getAssignedMaterialsIdentifiers() const {
    return this->pimpl->getAssignedMaterialsIdentifiers();
  }  // end of getAssignedMaterialsIdentifiers

  void NonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<AbstractBoundaryCondition> f) {
    this->pimpl->addBoundaryCondition(std::move(f));
  }  // end of addBoundaryCondition

  bool NonLinearEvolutionProblem::addBoundaryCondition(
      Context& ctx, std::unique_ptr<AbstractBoundaryCondition> f) noexcept {
    return this->pimpl->addBoundaryCondition(ctx, std::move(f));
  }  // end of addBoundaryCondition

  void NonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition> bc) {
    this->pimpl->addBoundaryCondition(std::move(bc));
  }  // end of NonLinearEvolutionProblem::addBoundaryCondition

  void NonLinearEvolutionProblem::addUniformDirichletBoundaryCondition(
      const Parameters& params) {
    this->addBoundaryCondition(
        std::make_unique<UniformDirichletBoundaryCondition>(*this, params));
  }  // end of addUniformBoundaryCondition

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

  std::map<size_type, size_type>
  NonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
                                                    const Parameter& m,
                                                    const std::string& l,
                                                    const std::string& b) {
    return this->pimpl->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  OptionalReference<const Material> NonLinearEvolutionProblem::getMaterial(
      Context& ctx, const Parameter& m, size_type b) const noexcept {
    return this->pimpl->getMaterial(ctx, m, b);
  }  // end of getMaterial

  OptionalReference<Material> NonLinearEvolutionProblem::getMaterial(
      Context& ctx, const Parameter& m, size_type b) noexcept {
    return this->pimpl->getMaterial(ctx, m, b);
  }  // end of getMaterial

  OptionalReference<const AbstractBehaviourIntegrator>
  NonLinearEvolutionProblem::getBehaviourIntegrator(
      Context& ctx, const Parameter& m, size_type b) const noexcept {
    return this->pimpl->getBehaviourIntegrator(ctx, m, b);
  }  // end of getBehaviourIntegrator

  OptionalReference<AbstractBehaviourIntegrator>
  NonLinearEvolutionProblem::getBehaviourIntegrator(Context& ctx,
                                                    const Parameter& m,
                                                    size_type b) noexcept {
    return this->pimpl->getBehaviourIntegrator(ctx, m, b);
  }  // end of getBehaviourIntegrator

  const Material& NonLinearEvolutionProblem::getMaterial(
      const Parameter& m) const {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  Material& NonLinearEvolutionProblem::getMaterial(const Parameter& m) {
    return this->pimpl->getMaterial(m);
  }  // end of getMaterial

  const AbstractBehaviourIntegrator&
  NonLinearEvolutionProblem::getBehaviourIntegrator(const size_type m) const {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  AbstractBehaviourIntegrator&
  NonLinearEvolutionProblem::getBehaviourIntegrator(const size_type m) {
    return this->pimpl->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblem::update() {
    this->pimpl->update();
  }  // end of update

  void NonLinearEvolutionProblem::revert() {
    this->pimpl->revert();
  }  // end of revert

  void NonLinearEvolutionProblem::setup(const real t, const real dt) {
    this->pimpl->setup(t, dt);
  }  // end of setup

  template <bool parallel>
  static NonLinearEvolutionProblemImplementation<parallel>&
  getImplementationInternal(AbstractNonLinearEvolutionProblem* const p) {
    auto* const pi =
        dynamic_cast<NonLinearEvolutionProblemImplementation<parallel>*>(p);
    if (pi == nullptr) {
      raise(
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
      raise(
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
    raise(
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
    raise(
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
      raise(
          "buildFacesDescription: "
          "unsupported parallel computations");
#endif
    }
    return buildFacesDescription(p.getImplementation<false>(), bid);
  }  // end of buildFacesDescription

  std::vector<std::pair<size_type, std::vector<std::vector<size_type>>>>
  getElementsDegreesOfFreedomOnBoundary(NonLinearEvolutionProblem& p,
                                        const size_type bid) {
    auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return getElementsDegreesOfFreedomOnBoundary(p.getImplementation<true>(),
                                                   bid);
#else
      raise(
          "getElementsDegreesOfFreedomOnBoundary: "
          "unsupported parallel computations");
#endif
    }
    return getElementsDegreesOfFreedomOnBoundary(p.getImplementation<false>(),
                                                 bid);
  }  // end of getElementsDegreesOfFreedomOnBoundary

  void computeResultantForceOnBoundary(
      mfem::Vector& F,
      NonLinearEvolutionProblem& p,
      const std::vector<
          std::pair<size_type, std::vector<std::vector<size_type>>>>&
          elts_dofs) {
    auto& fed = p.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      computeResultantForceOnBoundary(F, p.getImplementation<true>(),
                                      elts_dofs);
#else
      raise(
          "computeResultantForceOnBoundary: "
          "unsupported parallel computations");
#endif
    } else {
      computeResultantForceOnBoundary(F, p.getImplementation<false>(),
                                      elts_dofs);
    }
  }  // end of computeResultantForceOnBoundary

}  // end of namespace mfem_mgis
