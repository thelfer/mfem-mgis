/*!
 * \file   src/NonLinearEvolutionProblemImplementationBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include <iostream>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"
#include "MFEMMGIS/AbstractDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  const char* const NonLinearEvolutionProblemImplementationBase::
      UseMultiMaterialNonLinearIntegrator =
          "UseMultiMaterialNonLinearIntegrator";

  std::vector<std::string>
  NonLinearEvolutionProblemImplementationBase::getParametersList() {
    return {NonLinearEvolutionProblemImplementationBase::
                UseMultiMaterialNonLinearIntegrator};
  }  // end of getParametersList

  [[nodiscard]] static MultiMaterialNonLinearIntegrator*
  buildMultiMaterialNonLinearIntegrator(
      attributes::Throwing,
      std::shared_ptr<FiniteElementDiscretization> fed,
      const Hypothesis h,
      const Parameters& p) {
    const auto* const n = NonLinearEvolutionProblemImplementationBase::
        UseMultiMaterialNonLinearIntegrator;
    if (contains(p, n)) {
      if (!get<bool>(throwing, p, n)) {
        return nullptr;
      }
    }
    return new MultiMaterialNonLinearIntegrator(fed, h);
  }  // end of buildMultiMaterialNonLinearIntegrator

  NonLinearEvolutionProblemImplementationBase::
      NonLinearEvolutionProblemImplementationBase(
          std::shared_ptr<FiniteElementDiscretization> fed,
          const Hypothesis h,
          const Parameters& p)
      : fe_discretization(fed),
        u0(getTrueVSize(*fed)),
        u1(getTrueVSize(*fed)),
        mgis_integrator(
            buildMultiMaterialNonLinearIntegrator(throwing, fed, h, p)),
        hypothesis(h) {
    this->u0 = real{0};
    this->u1 = real{0};
  }  // end of NonLinearEvolutionProblemImplementationBase

  FiniteElementDiscretization& NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretization() noexcept {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  const FiniteElementDiscretization&
  NonLinearEvolutionProblemImplementationBase::getFiniteElementDiscretization()
      const noexcept {
    return *(this->fe_discretization);
  }  // end of getFiniteElementDiscretization

  std::shared_ptr<FiniteElementDiscretization>
  NonLinearEvolutionProblemImplementationBase::
      getFiniteElementDiscretizationPointer() noexcept {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  const mfem::Vector& NonLinearEvolutionProblemImplementationBase::
      getUnknownsAtBeginningOfTheTimeStep() const {
    return this->u0;
  }  // end of getUnknownsAtBeginningOfTheTimeStep

  mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep() {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  const mfem::Vector&
  NonLinearEvolutionProblemImplementationBase::getUnknownsAtEndOfTheTimeStep()
      const {
    return this->u1;
  }  // end of getUnknownsAtEndOfTheTimeStep

  mfem::Vector& NonLinearEvolutionProblemImplementationBase::getUnknowns(
      const TimeStepStage ts) noexcept {
    if (ts == bts) {
      return this->u0;
    }
    return this->u1;
  }  // end of getUnknowns

  const mfem::Vector& NonLinearEvolutionProblemImplementationBase::getUnknowns(
      const TimeStepStage ts) const noexcept {
    if (ts == bts) {
      return this->u0;
    }
    return this->u1;
  }  // end of getUnknowns

  void NonLinearEvolutionProblemImplementationBase::revert() {
    this->u1 = this->u0;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->revert();
    }
  }  // end of revert

  void NonLinearEvolutionProblemImplementationBase::update() {
    this->u0 = this->u1;
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->update();
    }
  }  // end of update

  [[nodiscard]] static bool checkMultiMaterialSupportEnabled(
      Context& ctx,
      const char* const n,
      const MultiMaterialNonLinearIntegrator* const p) {
    if (p == nullptr) {
      return ctx.registerErrorMessage(
          "NonLinearEvolutionProblemImplementationBase::" + std::string{n} +
          ": multi material support has been disabled");
    }
    return true;
  }  // end of checkMultiMaterialSupportEnabled

  static void checkMultiMaterialSupportEnabled(
      const char* const n, const MultiMaterialNonLinearIntegrator* const p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    checkMultiMaterialSupportEnabled(ctx, n, p) | or_raise;
  }  // end of checkMultiMaterialSupportEnabled

  bool NonLinearEvolutionProblemImplementationBase::setMaterialsNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->getFiniteElementDiscretization().setMaterialsNames(ctx, ids);
  }  // end of setMaterialsNames

  bool NonLinearEvolutionProblemImplementationBase::setBoundariesNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->getFiniteElementDiscretization().setBoundariesNames(ctx, ids);
  }  // end of setBoundariesNames

  void NonLinearEvolutionProblemImplementationBase::setMaterialsNames(
      const std::map<size_type, std::string>& ids) {
    this->getFiniteElementDiscretization().setMaterialsNames(ids);
  }

  void NonLinearEvolutionProblemImplementationBase::setBoundariesNames(
      const std::map<size_type, std::string>& ids) {
    this->getFiniteElementDiscretization().setBoundariesNames(ids);
  }

  std::optional<size_type>
  NonLinearEvolutionProblemImplementationBase::getMaterialIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    return this->getFiniteElementDiscretization().getMaterialIdentifier(ctx, p);
  }  // end of getMaterialIdentifier

  std::optional<size_type>
  NonLinearEvolutionProblemImplementationBase::getBoundaryIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    return this->getFiniteElementDiscretization().getBoundaryIdentifier(ctx, p);
  }  // end of getBoundaryIdentifier

  std::optional<std::vector<size_type>>
  NonLinearEvolutionProblemImplementationBase::getMaterialsIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    return this->getFiniteElementDiscretization().getMaterialsIdentifiers(ctx,
                                                                          p);
  }  // end of getMaterialsIdentifiers

  std::optional<std::vector<size_type>>
  NonLinearEvolutionProblemImplementationBase::getBoundariesIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    return this->getFiniteElementDiscretization().getBoundariesIdentifiers(ctx,
                                                                           p);
  }  // end of getBoundariesIdentifiers

  size_type NonLinearEvolutionProblemImplementationBase::getMaterialIdentifier(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getMaterialIdentifier(p);
  }  // end of getMaterialIdentifier

  size_type NonLinearEvolutionProblemImplementationBase::getBoundaryIdentifier(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getBoundaryIdentifier(p);
  }  // end of getBoundaryIdentifier

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getMaterialsIdentifiers(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getMaterialsIdentifiers(p);
  }  // end of getMaterialsIdentifiers

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getBoundariesIdentifiers(
      const Parameter& p) const {
    return this->getFiniteElementDiscretization().getBoundariesIdentifiers(p);
  }  // end of getBoundariesIdentifiers

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getAssignedMaterialsIdentifiers()
      const noexcept {
    if (this->mgis_integrator == nullptr) {
      return {};
    }
    return this->mgis_integrator->getAssignedMaterialsIdentifiers();
  }  // end of getAssignedMaterialsIdentifiers

  std::optional<std::map<size_type, size_type>>
  NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      Context& ctx,
      const std::string& n,
      const Parameter& m,
      const std::string& l,
      const std::string& b) noexcept {
    return this->addBehaviourIntegrator(ctx, n, m, l, b, {});
  }  // end of addBehaviourIntegrator

  std::optional<std::map<size_type, size_type>>
  NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      Context& ctx,
      const std::string& n,
      const Parameter& m,
      const std::string& l,
      const std::string& b,
      const Parameters& params) noexcept {
    checkMultiMaterialSupportEnabled("addBehaviourIntegrator",
                                     this->mgis_integrator);
    auto bids = std::map<size_type, size_type>{};
    for (const auto& mid : this->getMaterialsIdentifiers(m)) {
      const auto obid = this->mgis_integrator->addBehaviourIntegrator(
          ctx, n, mid, l, b, params);
      if (isInvalid(obid)) {
        return {};
      }
      bids.insert({mid, *obid});
    }
    return bids;
  }  // end of addBehaviourIntegrator

  std::map<size_type, size_type>
  NonLinearEvolutionProblemImplementationBase::addBehaviourIntegrator(
      const std::string& n,
      const Parameter& m,
      const std::string& l,
      const std::string& b) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return this->addBehaviourIntegrator(ctx, n, m, l, b) | or_raise;
  }  // end of addBehaviourIntegrator

  OptionalReference<const Material>
  NonLinearEvolutionProblemImplementationBase::getMaterial(
      Context& ctx, const Parameter& m, const size_type b) const noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "support for mgis integrator has been disabled");
    }
    const auto om = this->getMaterialIdentifier(ctx, m);
    if (isInvalid(om)) {
      return {};
    }
    return this->mgis_integrator->getMaterial(ctx, *om, b);
  }  // end of getMaterial

  OptionalReference<Material>
  NonLinearEvolutionProblemImplementationBase::getMaterial(
      Context& ctx, const Parameter& m, const size_type b) noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "support for mgis integrator has been disabled");
    }
    const auto om = this->getMaterialIdentifier(ctx, m);
    if (isInvalid(om)) {
      return {};
    }
    return this->mgis_integrator->getMaterial(ctx, *om, b);
  }  // end of getMaterial

  std::optional<size_type>
  NonLinearEvolutionProblemImplementationBase::getNumberOfBehaviourIntegrators(
      Context& ctx, const Parameter& m) const noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "support for mgis integrator has been disabled");
    }
    const auto om = this->getMaterialIdentifier(ctx, m);
    if (isInvalid(om)) {
      return {};
    }
    return this->mgis_integrator->getNumberOfBehaviourIntegrators(ctx, *om);
  }  // end of getNumberOfBehaviourIntegrators

  OptionalReference<const AbstractBehaviourIntegrator>
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      Context& ctx, const Parameter& m, const size_type b) const noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "support for mgis integrator has been disabled");
    }
    const auto om = this->getMaterialIdentifier(ctx, m);
    if (isInvalid(om)) {
      return {};
    }
    return this->mgis_integrator->getBehaviourIntegrator(ctx, *om, b);
  }  // end of getBehaviourIntegrator

  OptionalReference<AbstractBehaviourIntegrator>
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      Context& ctx, const Parameter& m, const size_type b) noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "support for mgis integrator has been disabled");
    }
    const auto om = this->getMaterialIdentifier(ctx, m);
    if (isInvalid(om)) {
      return {};
    }
    return this->mgis_integrator->getBehaviourIntegrator(ctx, *om, b);
  }  // end of getBehaviourIntegrator

  const Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const Parameter& m) const {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(this->getMaterialIdentifier(m));
  }  // end of getMaterial

  Material& NonLinearEvolutionProblemImplementationBase::getMaterial(
      const Parameter& m) {
    checkMultiMaterialSupportEnabled("getMaterial", this->mgis_integrator);
    return this->mgis_integrator->getMaterial(this->getMaterialIdentifier(m));
  }  // end of getMaterial

  const AbstractBehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) const {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  AbstractBehaviourIntegrator&
  NonLinearEvolutionProblemImplementationBase::getBehaviourIntegrator(
      const size_type m) {
    checkMultiMaterialSupportEnabled("getBehaviourIntegrator",
                                     this->mgis_integrator);
    return this->mgis_integrator->getBehaviourIntegrator(m);
  }  // end of getBehaviourIntegrator

  void NonLinearEvolutionProblemImplementationBase::setTimeIncrement(
      const real dt) {
    if (this->mgis_integrator != nullptr) {
      this->mgis_integrator->setTimeIncrement(dt);
    }
  }  // end of setTimeIncrement

  void NonLinearEvolutionProblemImplementationBase::setMacroscopicGradients(
      const std::vector<real>& g) {
    checkMultiMaterialSupportEnabled("setMacroscopicGradients",
                                     this->mgis_integrator);
    this->mgis_integrator->setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  bool NonLinearEvolutionProblemImplementationBase::setSolverParameters(
      Context& ctx, const Parameters& params) noexcept {
#ifdef MFEM_USE_PETSC
    if (usePETSc()) {
      return mfem_mgis::setSolverParameters(ctx, *(this->petsc_solver), params);
    }
    return mfem_mgis::setSolverParameters(ctx, *(this->solver), params);
#else  /* MFEM_USE_PETSC */
    return mfem_mgis::setSolverParameters(ctx, *(this->solver), params);
#endif /* MFEM_USE_PETSC */
  }    // end of setSolverParameters

  void NonLinearEvolutionProblemImplementationBase::setSolverParameters(
      const Parameters& params) {
#ifdef MFEM_USE_PETSC
    if (usePETSc()) {
      mfem_mgis::setSolverParameters(throwing, *(this->petsc_solver), params);
    } else {
      mfem_mgis::setSolverParameters(throwing, *(this->solver), params);
    }
#else  /* MFEM_USE_PETSC */
    mfem_mgis::setSolverParameters(throwing, *(this->solver), params);
#endif /* MFEM_USE_PETSC */
  }    // end of setSolverParameters

  std::vector<size_type>
  NonLinearEvolutionProblemImplementationBase::getEssentialDegreesOfFreedom()
      const {
    auto all_dofs = std::vector<mfem_mgis::size_type>{};
    for (const auto& bc : this->dirichlet_boundary_conditions) {
      auto dofs = bc->getHandledDegreesOfFreedom();
      all_dofs.insert(all_dofs.end(), dofs.begin(), dofs.end());
    }
    return all_dofs;
  }  // end of getEssentialDegreesOfFreedom

  [[nodiscard]] const std::vector<
      std::unique_ptr<AbstractDirichletBoundaryCondition>>&
  NonLinearEvolutionProblemImplementationBase::getDirichletBoundaryConditions()
      const noexcept {
    return this->dirichlet_boundary_conditions;
  }  // end of getDirichletBoundaryConditions

  const std::vector<std::unique_ptr<AbstractBoundaryCondition>>&
  NonLinearEvolutionProblemImplementationBase::getBoundaryConditions()
      const noexcept {
    return this->boundary_conditions;
  }  // end of getBoundaryConditions

  bool NonLinearEvolutionProblemImplementationBase::setup(
      Context& ctx, const real t, const real dt) noexcept {
    CatchTimeSection("NLEPIB::setup");
    const auto success = [this, &ctx, t, dt] {
      if (this->initialization_phase) {
        if (!this->dirichlet_boundary_conditions.empty()) {
          this->markDegreesOfFreedomHandledByDirichletBoundaryConditions(
              this->getEssentialDegreesOfFreedom());
        }
      }
      this->initialization_phase = false;
      for (const auto& bc : this->dirichlet_boundary_conditions) {
        bc->updateImposedValues(this->u1, t + dt);
      }
      for (const auto& bc : this->boundary_conditions) {
        bc->setup(t, dt);
      }
      if (this->mgis_integrator != nullptr) {
        if (!this->mgis_integrator->setup(ctx, t, dt)) {
          return false;
        }
      }
      return true;
    }();
    return isTrueOnAllProcesses(*(this->fe_discretization), success);
  }  // end of setup

  void NonLinearEvolutionProblemImplementationBase::setup(const real t,
                                                          const real dt) {
    auto ctx = Context{};
    if (!this->setup(ctx, t, dt)) {
      raise(ctx.getErrorMessage());
    }
  }  // end of setup

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s) {
    if (usePETSc()) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::updateLinearSolver: "
          "call to this method is meaningless if PETSc is used");
    }
    this->linear_solver_preconditioner.reset();
    this->linear_solver = std::move(s);
    this->solver->setLinearSolver(*(this->linear_solver));
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      std::unique_ptr<LinearSolver> s,
      std::unique_ptr<LinearSolverPreconditioner> p) {
    if (usePETSc()) {
      mgis::raise(
          "NonLinearEvolutionProblemImplementationBase::updateLinearSolver: "
          "call to this method is meaningless if PETSc is used");
    }
    if (p != nullptr) {
      auto* const isolver = dynamic_cast<IterativeSolver*>(s.get());
      if (isolver != nullptr) {
        isolver->SetPreconditioner(*p);
      }
      this->updateLinearSolver(std::move(s));
      this->linear_solver_preconditioner = std::move(p);
    } else {
      this->updateLinearSolver(std::move(s));
    }
  }  // end of updateLinearSolver

  void NonLinearEvolutionProblemImplementationBase::updateLinearSolver(
      LinearSolverHandler s) {
    CatchTimeSection("NLEPIB::updateLinearSolver");
    this->updateLinearSolver(std::move(s.linear_solver),
                             std::move(s.preconditioner));
  }  // end of updateLinearSolver

  std::optional<LinearizedOperators>
  NonLinearEvolutionProblemImplementationBase::getLinearizedOperators(
      Context& ctx, const mfem::Vector& U) noexcept {
    if (this->mgis_integrator == nullptr) {
      return ctx.registerErrorMessage(
          "no multiple material integrator defined");
    }
    return this->mgis_integrator->getLinearizedOperators(U);
  }  // end of getLinearizedOperators

  [[nodiscard]] bool NonLinearEvolutionProblemImplementationBase::
      areStiffnessOperatorsFromLastIterationAvailable() const noexcept {
    return this->hasStiffnessOperatorsBeenComputed;
  }  // end of areStiffnessOperatorsFromLastIterationAvailable

  void NonLinearEvolutionProblemImplementationBase::setPredictionPolicy(
      const PredictionPolicy& p) noexcept {
    this->prediction_policy = p;
  }  // end of setPredictionPolicy

  PredictionPolicy
  NonLinearEvolutionProblemImplementationBase::getPredictionPolicy()
      const noexcept {
    return this->prediction_policy;
  }  // end of getPredictionPolicy

  NonLinearResolutionOutput NonLinearEvolutionProblemImplementationBase::solve(
      const real t, const real dt) {
    auto ctx = Context{};
    const auto r = this->solve(ctx, t, dt);
    if (isInvalid(r)) {
      ctx.log() << ctx.getErrorMessage() << '\n';
    }
    return r;
  }  // end of solve

  NonLinearResolutionOutput NonLinearEvolutionProblemImplementationBase::solve(
      Context& ctx, const real t, const real dt) noexcept {
    CatchTimeSection("NLEPIB::solve");
    this->setTimeIncrement(dt);
    this->setup(t, dt);
    NonLinearResolutionOutput output;
    if (this->prediction_policy.strategy !=
        PredictionStrategy::DEFAULT_PREDICTION) {
      const auto onorm = this->computePrediction(ctx, t, dt);
      if (isInvalid(onorm)) {
        output.status = false;
        return output;
      }
      if (!usePETSc()) {
        if (!this->solver->setReferenceResidualNorm(ctx, *onorm)) {
          output.status = false;
          return output;
        }
      }
    }
    auto fill_output = [&output](auto& s) {
      output.status = s.GetConverged();
      output.iterations = s.GetNumIterations();
      output.final_residual_norm = s.GetFinalNorm();
    };
    if (usePETSc()) {
#ifdef MFEM_USE_PETSC
      mfem::Vector zero;
      // PETSc solver somehow "requires" zero as a first argument of Mult.
      // If it is not the case, one should take care of memory management.
      this->petsc_solver->Mult(zero, this->u1);
      fill_output(*(this->petsc_solver));
#else  /* MFEM_USE_PETSC */
      raise("Support for PETSc is deactivated");
#endif /* MFEM_USE_PETSC */
    } else {
      this->solver->setContext(ctx);
      this->solver->Mult(this->u0, this->u1);
      this->solver->unsetContext();
      fill_output(*(this->solver));
      output.initial_residual_norm = this->solver->GetInitialNorm();
      this->solver->unsetReferenceResidualNorm();
    }
    return output;
  }  // end of solve

  NonLinearEvolutionProblemImplementationBase::
      ~NonLinearEvolutionProblemImplementationBase() = default;

}  // end of namespace mfem_mgis
