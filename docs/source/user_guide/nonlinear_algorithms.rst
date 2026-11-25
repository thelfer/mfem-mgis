.. _mfem_mgis_nonlinear_algorithms:

==============================
Available nonlinear algorithms
==============================

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Prediction of the solution
==========================

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step, see below for details.

This can be changed by using the :cxx:`setPredictionPolicy`
method, as follows:

.. code-block:: c++

   // use the elastic operator by default
   mechanics.setPredictionPolicy(
      {.strategy = PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION});

Available strategies
--------------------

Default prediction (:cxx:`PredictionStrategy::DEFAULT_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step.

.. warning::

   In mechanics, this may lead to very high increments of the deformation
   gradients or the strain in the neighboring elements of boundaries where
   evolving displacements are imposed.
   
Prediction for the state at the beginning of the time step (:cxx:`PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The prediction for the state at the beginning of the time step strategy
determines the increment of the displacement :math:`\Delta\,\mathbb{u}`
by solving the following linear system:

.. math::

   \mathbb{K}\,\cdot\,\Delta\,\mathbb{u} = \ets{\mathbb{F}_{e}}-\bts{\mathbb{F}_{i}}

where:

- :math:`\mathbb{K}_{e}` denotes one of the prediction operator (see below).
- :math:`\bts{\mathbb{F}_{e}}` denotes the external forces at the beginning
  of the time step.
- :math:`\bts{\mathbb{F}_{i}}` denotes the inner forces at the beginning
  of the time step.
- :math:`\Delta\,\mathbb{u}` is submitted to the increment of the
  imposed Dirichlet boundary conditions.

.. note::

   Although the wording explicitly refers to mechanics, this equation
   applies to all physics.

The following prediction operators can be chosen:

- :cxx:`PredictionOperator::ELASTIC`: the elastic operator
- :cxx:`PredictionOperator::SECANT`: the secant operator is typically
  defined by the elastic operator
- :cxx:`PredictionOperator::TANGENT_PREDICTION`: the tangent operator,
  defined by the time-continuous derivative of the thermodynamic force
  with respect to the gradients.
- :cxx:`PredictionOperator::LAST_ITERATE_OPERATOR`: this operator reuses
  the one computed at the last iteration of the previous time step. At
  the first time step, the elastic operator is used.

.. note::

   Internally, the :cxx:`computePrediction` methods of the
   `NonLinearResolutionImplementation` classes computes the opposite of the 
   prediction.

Prediction based on a behaviour integration with constant gradients (:cxx:`PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cxx:`CONSTANT_GRADIENTS_INTEGRATION` strategy determines the
increment of the unknown :math:`\Delta\,{u}` by solving the following
linear system:
     
.. math::

   \tilde{\mathbb{K}}\,\cdot\,\Delta\,\mathbb{u} =
   \ets{\mathbb{F}_{e}}-\ets{\tilde{\mathbb{F}}_{i}}
     
where:
     
- :math:`\tilde{\mathbb{K}}` denotes the operator computed at the end of
  the behaviour integration.
- :math:`\bts{\mathbb{F}_{e}}` denotes the external forces at the
  beginning of the time step.
- :math:`\ets{\tilde{\mathbb{F}}_{i}}` denotes an approximation inner
  forces at the end of the time step computed by assuming that the
  gradients are constant over the time (and thus equal to their values at
  the beginning of the time step).
- :math:`\Delta\,\mathbb{u}` is submitted to the increment of the
  imposed Dirichlet boundary conditions.

The behaviour integration allows taking into account:

- the evolution of stress-free strain over the time step (thermal
  expansion, swelling, etc..),
- the viscoplastic relaxation of the stress. This relaxation can be
  discarded by integrating the behaviour with a null time step.

The following operators are available:

- :cxx:`IntegrationOperator::ELASTIC`: the elastic operator,
- :cxx:`IntegrationOperator::SECANT`: the secant operator is typically
  defined by the elastic operator affected by damage,
- :cxx:`IntegrationOperator::TANGENT`: the tangent operator, defined by
  the time-continuous derivative of the thermodynamic force with respect
  to the gradients,
- :cxx:`IntegrationOperator::CONSISTENT_TANGENT`: the consistent
  tangent operator, defined by the derivative of the thermodynamic force
  with respect to the gradients at the end of the time step. See
  :cite:`simo_consistent_1985` for details.

.. note::

   Internally, the :cxx:`computePrediction` methods of the
   `NonLinearResolutionImplementation` classes computes the opposite of the 
   prediction.

Newton algorithm
================

The default nonlinear solver is the Newton's algorithm.

Let :math:`\mathbb{u}^{(n)}` be the estimate of the solution at the
:math:`n^{\text{th}}` iteration. The Newton's algorithm determines a new
estimate :math:`\mathbb{u}^{(n+1)}` of the solution by computing a
correction :math:`\delta\,\mathbb{u}^{(n)}` as follows:

.. math::

   \mathbb{K}\,\cdot\,\delta\,\mathbb{u}^{(n)} = \ets{\mathbb{F}_{e}}^{(n)}-\ets{\mathbb{F}_{i}}^{(n)}

where:

- :math:`\mathbb{K}` denotes one of the search operator.
- :math:`\ets{\mathbb{F}_{e}}^{(n)}` denotes current estimate of the
  external forces at the end of the time step.
- :math:`\ets{\mathbb{F}_{i}}^{(n)}` denotes current estimate of the inner
  forces at the end of the time step.

Currently, the search operator is given by the consistent tangent operator.

The new estimate :math:`\mathbb{u}^{(n+1)}` is calculated as follows:

.. math::

   \mathbb{u}^{(n+1)} = \mathbb{u}^{(n)}+\delta\,\mathbb{u}^{(n)}

Handling of integration failure
-------------------------------

When the new estimate leads to an integration failure, a new estimate of
the solution is calculated by keeping the same search direction but by
reducing the amplitude of the correction by a factor :math:`2` as
follows:

.. math::

   \mathbb{u}^{(n+1)} = \mathbb{u}^{(n)}+\frac{1}{2}\,\delta\,\mathbb{u}^{(n)}

This reduction of the amplitude is repeated until an estimation of the
solution that does not lead to an integration failure is found.

.. _mfem_mgis_additional_convergence_criteria:

Additional convergence criteria
===============================

By default, the nonlinear solver (such as :cxx:`NewtonSolver`) considers a time
step converged when the residual norm satisfies:

.. math::

   \|\mathbb{r}\| \le \max\paren{\varepsilon_{r}\,\|\mathbb{r}_0\|,\; \varepsilon_{a}}

where :math:`\|\mathbb{r}_0\|` is the reference residual norm (obtained from the
prediction step if enabled, or at iteration ``0`` otherwise), and
:math:`\varepsilon_{r}` and :math:`\varepsilon_{a}` denote the relative and
absolute tolerances configured through :cxx:`setSolverParameters`.

In some simulations, this standard residual test alone is not sufficient. For
example, one may wish to:

- enforce an extra stopping condition before accepting a converged state;
- abort a resolution early when an undesirable state or divergence indicator
  is detected;
- dynamically adapt numerical, material, or behaviour parameters during the
  prediction step and across nonlinear iterations (continuation on material
  stiffness, adaptive linear solver tolerances, or damped initial tangent
  operators).

``MFEM/MGIS`` supports these use cases via the
:cxx:`AbstractAdditionalConvergenceCriterion` interface declared in
``MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx``.

Defining an additional convergence criterion
--------------------------------------------

To create a custom convergence criterion, derive a class from
:cxx:`mfem_mgis::AbstractAdditionalConvergenceCriterion` and implement its three
pure virtual methods:

- :cxx:`void helper() noexcept`:
  called during problem setup (:cxx:`NonLinearEvolutionProblemImplementationBase::setup`),
  **before** the prediction step is computed. It is typically used to ensure
  that parameters or material properties hold their nominal values for the
  prediction step.
- :cxx:`void reset() noexcept`:
  called during :cxx:`NonLinearEvolutionProblemImplementationBase::solve`,
  **after** the prediction step and right before the nonlinear solver runs. It
  is typically used to initialize or arm parameters for iteration ``0`` of the
  nonlinear algorithm.
- :cxx:`std::optional<bool> check(Context& ctx, const CheckArguments& args) noexcept`:
  called at every iteration of the nonlinear solver (starting at iteration
  ``0``), before checking the maximum iteration limit or computing the Newton
  correction.

The :cxx:`CheckArguments` structure passed to :cxx:`check` provides read-only
access to the current solver state:

- ``residual_norm`` (``const real``): current residual norm :math:`\|\mathbb{r}\|`.
- ``reference_residual_norm`` (``const real``): reference residual norm :math:`\|\mathbb{r}_0\|`.
- ``iter`` (``const int``): current iteration counter (starting at ``0``; also
  incremented by line-search amplitude-halving steps).
- ``max_iter`` (``const int``): maximum number of nonlinear iterations allowed.
- ``converged`` (``const bool``): result of the standard residual convergence
  test (``residual_norm <= norm_goal``).
- ``u`` (``const mfem::Vector&``): current estimate of the unknowns vector.

The return value of :cxx:`check` controls the solver loop as follows:

- returning ``true`` indicates that this additional criterion is satisfied;
- returning ``false`` **prevents** convergence at the current iteration even if
  ``args.converged`` is ``true``, causing the solver to perform another
  iteration (useful after updating continuation parameters);
- returning ``std::nullopt`` (for example via
  ``return ctx.registerErrorMessage("...");``) reports a failure, marks the
  nonlinear resolution as non-converged, and immediately breaks out of the
  solver loop.

.. note::

   An additional convergence criterion can only **strengthen** convergence or
   trigger an abort: it cannot force convergence when ``args.converged`` is
   ``false``. Specifically, the solver combines the standard test and all
   registered criteria using
   ``args.converged && crit_1 && crit_2 && ...``.

Registering criteria on a problem
---------------------------------

Criteria are registered on the problem's nonlinear solver through
:cxx:`getSolver()` and :cxx:`addAdditionalConvergenceCriterion`:

.. code-block:: c++

   #include "MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx"

   auto criterion = std::make_shared<MyConvergenceCriterion>(/* ... */);
   problem.getSolver().addAdditionalConvergenceCriterion(criterion);

Multiple criteria may be registered on the same solver. At each iteration,
:cxx:`check` is invoked on **every** registered criterion in registration order
without short-circuiting on ``false``, ensuring that parameter updates or
monitoring side effects occur for all criteria unless an earlier criterion
returns ``std::nullopt``.

Usage examples
--------------

Minimal custom check and error reporting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following example enforces a minimum number of iterations and demonstrates
reporting an error via :cxx:`Context`:

.. code-block:: c++

   #include "MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx"

   struct MinimumIterationsCriterion final
       : mfem_mgis::AbstractAdditionalConvergenceCriterion {
     explicit MinimumIterationsCriterion(const int min_it) noexcept
         : min_iterations(min_it) {}

     void helper() noexcept override {
       // Nothing to restore before prediction
     }

     void reset() noexcept override {
       this->checks_done = 0;
     }

     [[nodiscard]] std::optional<bool> check(
         mfem_mgis::Context& ctx,
         const CheckArguments& args) noexcept override {
       if (args.residual_norm > 1e12) {
         return ctx.registerErrorMessage(
             "residual norm diverged beyond admissible threshold");
       }
       ++(this->checks_done);
       if (!args.converged) {
         return false;
       }
       return args.iter >= this->min_iterations;
     }

     ~MinimumIterationsCriterion() noexcept override = default;

    private:
     int min_iterations;
     int checks_done = 0;
   };

Adaptive linear solver tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For large nonlinear problems solved with Krylov iterative linear solvers, solving
early Newton systems to high linear accuracy can be wasteful. A continuation
strategy starts each time step with a relaxed linear tolerance (e.g.,
:math:`100\,\text{tol}_0` in :cxx:`reset()`) and tightens it progressively once
the nonlinear residual satisfies the tolerance, until the target linear
tolerance :math:`\text{tol}_0` is reached:

.. code-block:: c++

   struct LinearToleranceCriterion final
       : mfem_mgis::AbstractAdditionalConvergenceCriterion {
     LinearToleranceCriterion(mfem_mgis::NonLinearEvolutionProblem& p,
                              mfem_mgis::Context& c,
                              const mfem_mgis::real target_tol,
                              mfem_mgis::Parameters ls_params,
                              std::string solver_name)
         : problem(p),
           ctx(c),
           tol0(target_tol),
           current_tol(target_tol),
           params(std::move(ls_params)),
           linear_solver_name(std::move(solver_name)) {}

     [[nodiscard]] bool updateLinearSolver(mfem_mgis::Context& c) noexcept {
       // Parameter name depends on the linear solver (e.g., "RelativeTolerance"
       // for MFEM Krylov solvers such as CGSolver or GMRESSolver)
       this->params.replaceOrInsert("RelativeTolerance", this->current_tol);
       return this->problem.setLinearSolver(
           c, this->linear_solver_name, this->params);
     }

     void helper() noexcept override {
       this->current_tol = this->tol0;
       this->update_failed = !this->updateLinearSolver(this->ctx);
     }

     void reset() noexcept override {
       this->current_tol = mfem_mgis::real{100} * this->tol0;
       this->update_failed = !this->updateLinearSolver(this->ctx);
     }

     [[nodiscard]] std::optional<bool> check(
         mfem_mgis::Context& c, const CheckArguments& args) noexcept override {
       if (this->update_failed) {
         return c.registerErrorMessage(
             "failed to update linear solver during setup/reset");
       }
       if (!args.converged) {
         return false;
       }
       if (this->current_tol <= this->tol0) {
         return true;
       }
       constexpr auto divisor = mfem_mgis::real{10};
       this->current_tol = std::max(this->tol0, this->current_tol / divisor);
       if (!this->updateLinearSolver(c)) {
         return {};
       }
       return false;
     }

     ~LinearToleranceCriterion() noexcept override = default;

    private:
     mfem_mgis::NonLinearEvolutionProblem& problem;
     mfem_mgis::Context& ctx;
     mfem_mgis::real tol0;
     mfem_mgis::real current_tol;
     mfem_mgis::Parameters params;
     std::string linear_solver_name;
     bool update_failed = false;
   };

Behaviour parameter damping and continuation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In stiff nonlinear problems such as third-medium contact, a criterion can also
adjust constitutive behaviour parameters dynamically during Newton iterations
using :cxx:`mgis::behaviour::setParameter`:

.. code-block:: c++

   struct TangentDampingCriterion final
       : mfem_mgis::AbstractAdditionalConvergenceCriterion {
     TangentDampingCriterion(mfem_mgis::NonLinearEvolutionProblem& p,
                             mfem_mgis::Context& c,
                             const mfem_mgis::real a_max,
                             std::vector<std::string> mats)
         : problem(p),
           ctx(c),
           alpha_max(a_max),
           materials(std::move(mats)) {}

     [[nodiscard]] bool setAlpha(mfem_mgis::Context& c,
                                 const mfem_mgis::real val) noexcept {
       for (const auto& name : this->materials) {
         auto om = this->problem.getMaterial(c, name, 0);
         if (mfem_mgis::isInvalid(om)) {
           return false;
         }
         mgis::behaviour::setParameter(om->b, "alpha_tangent", val);
       }
       return true;
     }

     void helper() noexcept override {
       // Keep nominal tangent (alpha = 0) during prediction
       this->update_failed = !this->setAlpha(this->ctx, 0);
     }

     void reset() noexcept override {
       // Arm damping for iteration 0 only
       this->prev_iter = 0;
       this->update_failed = !this->setAlpha(this->ctx, this->alpha_max);
     }

     [[nodiscard]] std::optional<bool> check(
         mfem_mgis::Context& c, const CheckArguments& args) noexcept override {
       if (this->update_failed) {
         return c.registerErrorMessage(
             "failed to update behaviour parameter during setup/reset");
       }
       if (args.converged) {
         this->prev_iter = 0;
         if (!this->setAlpha(c, 0)) {
           return {};
         }
         return true;
       }
       // Re-arm damping for one step if line-search halving occurred
       const bool had_halving =
           (args.iter > 0) && (args.iter > this->prev_iter + 1);
       this->prev_iter = args.iter;
       if (!this->setAlpha(c, had_halving ? this->alpha_max : mfem_mgis::real{0})) {
         return {};
       }
       return false;
     }

     ~TangentDampingCriterion() noexcept override = default;

    private:
     mfem_mgis::NonLinearEvolutionProblem& problem;
     mfem_mgis::Context& ctx;
     mfem_mgis::real alpha_max;
     int prev_iter = 0;
     std::vector<std::string> materials;
     bool update_failed = false;
   };

Important execution notes
-------------------------

1. **One-iteration lag for residual and Jacobian updates**:
   At the beginning of :cxx:`NewtonSolver::Mult` (before ``iter == 0``) and at
   the end of each iteration (after updating ``x``), the solver integrates the
   constitutive behaviour and computes the residual before calling :cxx:`check`.
   Consequently:

   - changes applied in :cxx:`helper()` affect the prediction step;
   - changes applied in :cxx:`reset()` affect iteration ``0``;
   - updates to constitutive behaviour parameters or external fields inside
     :cxx:`check()` affect the **next** integration (iteration ``iter + 1``);
   - updates to the linear solver via :cxx:`problem.setLinearSolver` inside
     :cxx:`check()` take effect immediately in :cxx:`computeNewtonCorrection` of
     the **current** iteration ``iter``.

2. **Uniform scalar material properties vs. fields and behaviour parameters**:
   During :cxx:`NonLinearEvolutionProblemImplementationBase::setup` (which calls
   :cxx:`BehaviourIntegratorBase::setup` *before* :cxx:`helper()`), uniform
   scalar material properties set with
   ``mgis::behaviour::setMaterialProperty(ctx, m.s1, name, scalar_value)`` are
   copied by value into the integrator's internal workspace. Modifying a uniform
   scalar material property inside :cxx:`helper()`, :cxx:`reset()`, or
   :cxx:`check()` therefore takes effect on the next call to ``setup`` (the next
   time step). For immediate within-step updates, prefer behaviour parameters
   (:cxx:`mgis::behaviour::setParameter`) or spatially-indexed buffers/evaluators
   (``std::span<real>`` or :cxx:`QPEvaluator`) bound to ``m.s1``.

3. **PETSc nonlinear solvers**:
   Additional convergence criteria rely on :cxx:`AbstractNonLinearSolver`
   (such as :cxx:`NewtonSolver`) and are not invoked when PETSc nonlinear
   solvers are enabled (``--use-petsc``).

4. **Internal architecture**:
   For a detailed description of how these hooks are wired into the nonlinear
   problem and solver classes, see
   :ref:`mfem_mgis_developer_guide_nonlinear_solvers`.