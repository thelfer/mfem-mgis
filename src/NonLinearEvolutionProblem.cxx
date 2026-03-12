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
#include "MFEMMGIS/AbstractDirichletBoundaryCondition.hxx"
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
    const auto h = mgis::behaviour::fromString(get<std::string>(
        throwing, p, NonLinearEvolutionProblem::HypothesisParameter));
    const auto fed = std::make_shared<FiniteElementDiscretization>(
        extract(throwing, p, FiniteElementDiscretization::getParametersList()));
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

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(MeshDiscretization& m,
                                                       const Hypothesis h,
                                                       const Parameters& p)
      : NonLinearEvolutionProblem(
            std::make_shared<FiniteElementDiscretization>(
                m,
                extract(throwing,
                        p,
                        FiniteElementDiscretization::getParametersList())),
            h,
            p) {}  // end of NonLinearEvolutionProblem

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(MeshDiscretization& m,
                                                       const Parameters& p)
      : NonLinearEvolutionProblem(
            std::make_shared<FiniteElementDiscretization>(
                m,
                extract(throwing,
                        p,
                        FiniteElementDiscretization::getParametersList())),
            mgis::behaviour::fromString(get<std::string>(
                throwing, p, NonLinearEvolutionProblem::HypothesisParameter)),
            remove(p, {NonLinearEvolutionProblem::HypothesisParameter})) {
  }  // end of NonLinearEvolutionProblem

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Parameters& p)
      : NonLinearEvolutionProblem(
            fed,
            mgis::behaviour::fromString(get<std::string>(
                throwing, p, NonLinearEvolutionProblem::HypothesisParameter)),
            remove(p, {NonLinearEvolutionProblem::HypothesisParameter})) {
  }  // end of NonLinearEvolutionProblem

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    checkParameters(throwing, p,
                    NonLinearEvolutionProblem::getParametersList());
    if (contains(p, NonLinearEvolutionProblem::HypothesisParameter)) {
      raise("parameter '" +
            std::string{NonLinearEvolutionProblem::HypothesisParameter} +
            "' is not allowed");
    }
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
  NonLinearEvolutionProblem::getFiniteElementDiscretization() noexcept {
    return this->pimpl->getFiniteElementDiscretization();
  }  // end of getFiniteElementDiscretization

  const FiniteElementDiscretization&
  NonLinearEvolutionProblem::getFiniteElementDiscretization() const noexcept {
    return this->pimpl->getFiniteElementDiscretization();
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblem::getFiniteElementDiscretizationPointer() noexcept {
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

  bool NonLinearEvolutionProblem::setSolverParameters(
      Context& ctx, const Parameters& params) noexcept {
    return this->pimpl->setSolverParameters(ctx, params);
  }  // end of setSolverParameters

  void NonLinearEvolutionProblem::setSolverParameters(
      const Parameters& params) {
    this->pimpl->setSolverParameters(params);
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

  bool NonLinearEvolutionProblem::setLinearSolver(
      Context& ctx, std::string_view n, const Parameters& params) noexcept {
    return this->pimpl->setLinearSolver(ctx, n, params);
  }  // end of setLinearSolver

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

  const std::vector<std::unique_ptr<AbstractDirichletBoundaryCondition>>&
  NonLinearEvolutionProblem::getDirichletBoundaryConditions() const noexcept {
    return this->pimpl->getDirichletBoundaryConditions();
  }  // end of getDirichletBoundaryConditions

  const std::vector<std::unique_ptr<AbstractBoundaryCondition>>&
  NonLinearEvolutionProblem::getBoundaryConditions() const noexcept {
    return this->pimpl->getBoundaryConditions();
  }  // end of getBoundaryConditions

  NonLinearResolutionOutput NonLinearEvolutionProblem::solve(
      Context& ctx, const real t, const real dt) noexcept {
    CatchTimeSection("NLEP::solve");
    if (!this->setup(ctx, t, dt)) {
      return InvalidResult{};
    }
    return this->pimpl->solve(ctx, t, dt);
  }  // end of solve

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
  NonLinearEvolutionProblem::getAssignedMaterialsIdentifiers() const noexcept {
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

  bool NonLinearEvolutionProblem::addBoundaryCondition(
      Context& ctx,
      std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept {
    return this->pimpl->addBoundaryCondition(ctx, std::move(bc));
  }  // end of NonLinearEvolutionProblem::addBoundaryCondition

  void NonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<AbstractDirichletBoundaryCondition> bc) {
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

  bool NonLinearEvolutionProblem::addPostProcessing(
      Context& ctx, std::string_view n, const Parameters& p) noexcept {
    return this->pimpl->addPostProcessing(ctx, n, p);
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

  std::optional<size_type>
  NonLinearEvolutionProblem::getNumberOfBehaviourIntegrators(
      Context& ctx, const Parameter& m) const noexcept {
    return this->pimpl->getNumberOfBehaviourIntegrators(ctx, m);
  }  // end of getNumberOfBehaviourIntegrators

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

  bool NonLinearEvolutionProblem::setup(Context& ctx,
                                        const real t,
                                        const real dt) noexcept {
    return this->pimpl->setup(ctx, t, dt);
  }  // end of setup

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
      reportUnsupportedParallelComputations();
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
      reportUnsupportedParallelComputations();
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
      reportUnsupportedParallelComputations();
#endif
    } else {
      computeResultantForceOnBoundary(F, p.getImplementation<false>(),
                                      elts_dofs);
    }
  }  // end of computeResultantForceOnBoundary

  [[nodiscard]] static bool resolveMaterialDependencies(
      Context& ctx,
      AbstractBehaviourIntegrator& bi,
      const AbstractBehaviourIntegrator& provider) {
    using mgis::behaviour::MaterialStateManager;
    if ((!bi.hasMaterial()) || (!provider.hasMaterial())) {
      return true;
    }
    const auto om = bi.getMaterial(ctx);
    const auto oprovider_material = provider.getMaterial(ctx);
    if (isInvalid(om) || isInvalid(oprovider_material)) {
      return false;
    }
    // material properties
    for (const auto& mp : om->b.mps) {
      if (!contains(oprovider_material->b.isvs, mp.name)) {
        continue;
      }
      const auto oiv = getVariable(ctx, oprovider_material->b.isvs, mp.name);
      if (isInvalid(oiv)) {
        return false;
      }
      if (mp.type_identifier != oiv->type_identifier) {
        return ctx.registerErrorMessage("incompatible type for variable '" +
                                        mp.name + "'");
      }
      try {
        const auto ovalues_bts =
            getInternalStateVariable(ctx, *oprovider_material, mp.name,
                                     Material::BEGINNING_OF_TIME_STEP);
        const auto ovalues_ets = getInternalStateVariable(
            ctx, *oprovider_material, mp.name, Material::END_OF_TIME_STEP);
        //         setMaterialProperty(m.s0, mp.name, values_bts.getValues(),
        //                             MaterialStateManager::EXTERNAL_STORAGE,
        //                             MaterialStateManager::NOUPDATE);
        //         setMaterialProperty(m.s1, mp.name, values_ets.getValues(),
        //                             MaterialStateManager::EXTERNAL_STORAGE,
        //                             MaterialStateManager::NOUPDATE);
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx);
      }
    }
    // external state variables
    for (const auto& esv : om->b.esvs) {
    }
    return true;
  }  // end of resolveBehaviourIntegratorsDependencies

  [[nodiscard]] static bool resolveBehaviourIntegratorsDependencies(
      Context& ctx,
      NonLinearEvolutionProblemImplementationBase& p,
      const NonLinearEvolutionProblemImplementationBase& provider) noexcept {
    if (&p == &provider) {
      return true;
    }
    const auto provider_materials = provider.getAssignedMaterialsIdentifiers();
    for (const auto& m : p.getAssignedMaterialsIdentifiers()) {
      if (std::find(provider_materials.begin(), provider_materials.end(), m) ==
          provider_materials.end()) {
        continue;
      }
      const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, m);
      const auto onbis2 = provider.getNumberOfBehaviourIntegrators(ctx, m);
      if (!areValid(onbis, onbis2)) {
        return false;
      }
      for (auto n = size_type{}; n != *onbis; ++n) {
        auto obi = p.getBehaviourIntegrator(ctx, m, n);
        if (isInvalid(obi)) {
          return false;
        }
        for (auto n2 = size_type{}; n2 != *onbis; ++n2) {
          const auto obi2 = provider.getBehaviourIntegrator(ctx, m, n2);
          if (isInvalid(obi2)) {
            return false;
          }
          if (!resolveMaterialDependencies(ctx, *obi, *obi2)) {
            return false;
          }
        }
      }
    }
    return true;
  }  // end of resolveBehaviourIntegratorsDependencies

  bool resolveBehaviourIntegratorsDependencies(
      Context& ctx,
      NonLinearEvolutionProblem& p,
      const NonLinearEvolutionProblem& provider) noexcept {
    const auto& fed = p.getFiniteElementDiscretization();
    if (static_cast<const MeshDiscretization&>(fed) !=
        static_cast<const MeshDiscretization&>(
            provider.getFiniteElementDiscretization())) {
      return ctx.registerErrorMessage("inconsistent meshes");
    }
    if (fed.describesAParallelComputation()) {
      // since meshes are the same, both problems describes a parallel
      // computation
#ifdef MFEM_USE_MPI
      return resolveBehaviourIntegratorsDependencies(
          ctx, p.getImplementation<true>(), provider.getImplementation<true>());
#else
      reportUnsupportedParallelComputations();
#endif
    }
    return resolveBehaviourIntegratorsDependencies(
        ctx, p.getImplementation<false>(), provider.getImplementation<false>());
  }  // end of resolveBehaviourIntegratorsDependencies

}  // end of namespace mfem_mgis
