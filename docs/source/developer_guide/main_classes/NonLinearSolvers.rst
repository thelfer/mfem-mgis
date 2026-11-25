.. _mfem_mgis_developer_guide_nonlinear_solvers:

===================================================
Nonlinear solvers and additional convergence checks
===================================================

.. contents:: Table of Contents
   :depth: 3
   :local:

This page describes the internal architecture of the nonlinear solver layer in
``MFEM/MGIS`` and explains how additional convergence criteria
(:cxx:`AbstractAdditionalConvergenceCriterion`) are integrated across the
evolution problem and solver classes.

Architectural overview
======================

To decouple problem management from concrete nonlinear algorithms (and allow
solvers beyond :cxx:`NewtonSolver` to be registered via
:cxx:`NonLinearSolverFactory`), the nonlinear solver interface and additional
convergence check dispatching are placed on :cxx:`AbstractNonLinearSolver` and
:cxx:`NonLinearSolverBase`, rather than solely on :cxx:`NewtonSolver`.

The responsibilities are distributed across the following classes:

1. :cxx:`AbstractAdditionalConvergenceCriterion`
   (``include/MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx``):
   defines the :cxx:`CheckArguments` payload and the pure virtual hooks
   :cxx:`helper()`, :cxx:`reset()`, and :cxx:`check(Context&, const CheckArguments&)`.
2. :cxx:`AbstractNonLinearSolver`
   (``include/MFEMMGIS/NonLinearSolvers/AbstractNonLinearSolver.hxx``,
   ``src/AbstractNonLinearSolver.cxx``):
   inherits from ``mfem::IterativeSolver`` and declares the polymorphic
   interface for registering criteria (:cxx:`addAdditionalConvergenceCriterion`)
   and dispatching the lifecycle hooks
   (:cxx:`processAdditionalConvergenceCriterionHelper`,
   :cxx:`processAdditionalConvergenceCriterionReset`,
   :cxx:`processAdditionalConvergenceCriterionCheck`).
3. :cxx:`NonLinearSolverBase`
   (``include/MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx``,
   ``src/NonLinearSolverBase.cxx``):
   provides the shared storage and ``final`` dispatch implementations for all
   native nonlinear solvers:

   - stores registered criteria in ``acc_actions``
     (``std::vector<std::shared_ptr<AbstractAdditionalConvergenceCriterion>>``);
   - implements :cxx:`addAdditionalConvergenceCriterion` (profiled under
     ``"NS::addAdditionalConvergenceCriterion"``);
   - implements :cxx:`processAdditionalConvergenceCriterionHelper` and
     :cxx:`processAdditionalConvergenceCriterionReset` by forwarding to each
     criterion in registration order;
   - implements :cxx:`processAdditionalConvergenceCriterionCheck` (profiled under
     ``"NS::processAdditionalConvergenceCriterionCheck"``);
   - registers the default unknowns-estimate action (calling
     ``p.integrate(u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR, {})``)
     in its constructor.
4. :cxx:`NewtonSolver`
   (``include/MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx``,
   ``src/NewtonSolver.cxx``):
   derives from :cxx:`NonLinearSolverBase` and invokes
   :cxx:`processAdditionalConvergenceCriterionCheck` inside
   :cxx:`NewtonSolver::Mult` at every Newton iteration.
5. :cxx:`NonLinearSolverFactory`
   (``include/MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx``,
   ``src/NonLinearSolverFactory.cxx``):
   registers built-in nonlinear solvers (currently ``"Newton"`` mapped to
   :cxx:`NewtonSolver` via :cxx:`StandardNonLinearSolverGenerator`).
6. :cxx:`AbstractNonLinearEvolutionProblem` / :cxx:`NonLinearEvolutionProblem` /
   :cxx:`NonLinearEvolutionProblemImplementationBase` /
   :cxx:`NonLinearEvolutionProblemImplementation<parallel>`:
   expose :cxx:`getSolver()` returning :cxx:`AbstractNonLinearSolver&`,
   instantiate ``this->solver`` (via :cxx:`NonLinearSolverFactory` when the
   ``NonLinearSolver`` parameter is provided, or :cxx:`NewtonSolver` by
   default), and trigger the ``Helper`` and ``Reset`` hooks during ``setup``
   and ``solve``.

.. note::

   When ``usePETSc()`` is active, :cxx:`NonLinearEvolutionProblemImplementation`
   instantiates ``this->petsc_solver`` (``mfem::PetscNonlinearSolver``) instead
   of ``this->solver``. In :cxx:`NonLinearEvolutionProblemImplementationBase::setup`
   and :cxx:`solve`, the additional convergence criterion hooks are guarded by
   ``if (this->solver != nullptr)``.

Evaluation semantics in :cxx:`NonLinearSolverBase`
==================================================

In ``src/NonLinearSolverBase.cxx``,
:cxx:`NonLinearSolverBase::processAdditionalConvergenceCriterionCheck` evaluates
all registered criteria using the following rules:

.. code-block:: c++

   bool cv = s.converged;
   for (auto& a : this->acc_actions) {
     std::optional<bool> result = a->check(ctx, s);
     if (isInvalid(result)) {
       return {};
     }
     cv = cv && *result;
   }
   return cv;

Key properties of this loop:

- **No short-circuiting on** ``false``: even when ``s.converged`` is
  ``false`` or an earlier criterion returned ``false``, subsequent criteria in
  ``acc_actions`` still have their :cxx:`check` method called so that any
  parameter updates or side effects take place consistently.
- **Immediate abort on invalid optional**: if any criterion returns
  ``std::nullopt`` (:cxx:`InvalidResult`), the loop exits immediately and
  returns an empty ``std::optional<bool>``, signaling an error to the solver.
- **Standard residual check cannot be bypassed**: ``cv`` starts from
  ``s.converged`` (``norm <= norm_goal``), so a criterion can prevent
  convergence (by returning ``false``) or abort (by returning ``std::nullopt``),
  but cannot force convergence while ``s.converged`` is ``false``.

Chronological call flow ("what is called when")
===============================================

The full call sequence across a :cxx:`Simulation` or direct calls to
:cxx:`NonLinearEvolutionProblem::solve` is summarized below:

.. code-block:: text

   [Initialization (once)]
     NonLinearEvolutionProblem(...)
       `-> NonLinearEvolutionProblemImplementation<parallel>(...)
             `-> instantiates this->solver (NewtonSolver or NonLinearSolverFactory)
     problem.getSolver().addAdditionalConvergenceCriterion(crit)
       `-> NonLinearSolverBase::addAdditionalConvergenceCriterion(crit)
             `-> stores crit in acc_actions

   [Per time step: problem.solve(ctx, t, dt)]
     1. NonLinearEvolutionProblem::setup(ctx, t, dt)
          `-> NonLinearEvolutionProblemImplementationBase::setup(ctx, t, dt)
                |-- setTimeIncrement(dt)
                |-- mark essential Dirichlet DOFs (first call only)
                |-- bc->updateImposedValues(u1, t + dt) and bc->setup(ctx, t, dt)
                |-- mgis_integrator->setup(ctx, t, dt)
                `-- solver->processAdditionalConvergenceCriterionHelper()
                      `-> for each a in acc_actions: a->helper()
          (PeriodicNonLinearEvolutionProblem::setup then calls setMacroscopicGradients)

     2. NonLinearEvolutionProblemImplementationBase::solve(ctx, t, dt)
          |-- if (prediction_policy.strategy != DEFAULT_PREDICTION):
          |     computePrediction(ctx, t, dt)   // uses state after helper()
          |     solver->setReferenceResidualNorm(ctx, *onorm)
          |-- solver->processAdditionalConvergenceCriterionReset()
          |     `-> for each a in acc_actions: a->reset()
          |-- solver->setContext(ctx)
          |-- solver->Mult(u0, u1)              // NewtonSolver::Mult
          |-- solver->unsetContext()
          `-- fill_output & solver->unsetReferenceResidualNorm()

   [Inside NewtonSolver::Mult(b, x)]
     Pre-loop (using state active after reset()):
       |-- processNewUnknownsEstimate(x)       // integrates behaviour (consistent tangent)
       |-- norm = updateResidual()             // evaluates residual r and ||r||
       |-- sets reference_residual_norm if unset (DEFAULT_PREDICTION)
       `-- norm_goal = max(rel_tol * ||r_0||, abs_tol), it = 0

     while (true):
       |-- Monitor(it, norm, r, x)
       |-- result = processAdditionalConvergenceCriterionCheck(ctx, CheckArguments{...})
       |     `-> calls a->check(ctx, args) for each criterion in acc_actions
       |-- if (isInvalid(result)) -> converged = 0; break;
       |-- converged = *result ? 1 : 0; if (converged) break;
       |-- if (it >= max_iter) break;
       |-- computeNewtonCorrection(c, r, x)    // assembles Jacobian, solves linear system
       |-- x -= c
       |-- processNewUnknownsEstimate(x)       // integrates behaviour (or halves c on failure)
       |-- updateResidual()                    // updates residual and norm
       `-- record LinearSolver info, ++it

Execution timing details
========================

Prediction vs. iteration ``0`` (:cxx:`helper` vs. :cxx:`reset`)
---------------------------------------------------------------

- :cxx:`helper()` executes inside ``setup`` **before** :cxx:`computePrediction`,
  allowing criteria to restore nominal parameters (e.g., nominal stiffness or
  nominal linear solver tolerance) for the prediction linear solve.
- :cxx:`reset()` executes inside ``solve`` **after** :cxx:`computePrediction`
  and **before** :cxx:`NewtonSolver::Mult` calls
  :cxx:`processNewUnknownsEstimate(x)`. Therefore, parameters configured in
  :cxx:`reset()` govern the initial behaviour integration, initial residual
  computation, and first Jacobian assembly at ``iter == 0``.

Inside the Newton loop (:cxx:`check`)
-------------------------------------

- :cxx:`check` is called **before** ``it >= max_iter`` and **before**
  :cxx:`computeNewtonCorrection`.
  - When a prediction already satisfies ``norm <= norm_goal`` at ``iter == 0``,
    :cxx:`check` is invoked once with ``args.iter == 0`` and
    ``args.converged == true``.
  - Even at ``iter == max_iter``, :cxx:`check` is still invoked before the
    maximum-iteration exit test.
- **Immediate effect vs. one-iteration lag**:
  - Calling :cxx:`problem.setLinearSolver` inside :cxx:`check` updates
    ``this->prec`` before :cxx:`computeNewtonCorrection(c, r, x)` runs on the
    **same** iteration ``it``.
  - Updating constitutive behaviour parameters (e.g., via
    :cxx:`mgis::behaviour::setParameter`) or field data inside :cxx:`check`
    does **not** re-run :cxx:`processNewUnknownsEstimate` for iteration ``it``;
    instead, it takes effect during :cxx:`processNewUnknownsEstimate(x)` at the
    end of iteration ``it``, affecting the residual and Jacobian evaluated at
    iteration ``it + 1``.
- **Interaction with** :cxx:`Simulation` **sub-stepping**:
  If :cxx:`check` returns ``std::nullopt`` (or convergence fails after
  ``max_iter``), :cxx:`NonLinearEvolutionProblem::solve` returns a
  :cxx:`NonLinearResolutionOutput` with ``status == false``. Inside
  :cxx:`Simulation::simulateOverATimeStep`, this maps to
  ``ExitStatus::recoverableError``, causing :cxx:`Simulation` to invoke
  ``problem.revert(ctx)``, reduce the time increment, and restart the step
  through ``setup`` (:cxx:`helper()`) and ``solve`` (:cxx:`reset()`).
