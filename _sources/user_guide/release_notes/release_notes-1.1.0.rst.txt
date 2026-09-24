.. _mfem_mgis_release_notes_1_1:

============================
Release notes of Version 1.1
============================

This version inherits from all the features introduced in:

- Version 1.0.1
- Version 1.0.2
- Version 1.0.3
- :ref:`mfem_mgis_release_notes_1_0_4`

.. contents:: Table of Contents
   :depth: 3
   :local:
   :backlinks: none

Highlights
==========

- `NonLinearEvolutionProblem` are now able to make a prediction of the
  solution, removing major convergence issues when imposed displacements
  are imposed.
- The :math:`\bar{F}` method, or FBar formulation, has been implemented
  for finite strain behaviours to handle nearly incompressible materials.
- The regularization proposed by Faltus et al. in the context of the
  third medium contact has been implemented for plane strain, plane
  stress and tridmensional hypotheses.
- `MGIS`'s contexts now handle gathering computation time information to
  create the performance table instead of the previously used
  `CatchTimeSection`.
- Many methods have been deprecated to have a consistent error handling
  scheme based on `MGIS`'s one. As such, many methods and functions now
  takes and `MGIS`'s :cxx:`Context` as their first argument.
- The name of the materials and boundaries are automatically retrieved
  from |MFEM|'s mesh.

Known incompatibilites
======================

- In previous versions, the failure of Hypre's linear iterative solvers
  were discarded by `MFEM/MGIS`'s Newton solver due to the lack of
  methods to test their convergence in |MFEM|'s version prior to 4.10.
  The parameter `DiscardLinearSolverFailure` can be passed to
  `MFEM/MGIS`'s Newton solver to recover the behavior of previous
  versions.

New features
============

Prediction of the solution
--------------------------

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
^^^^^^^^^^^^^^^^^^^^

Default prediction (:cxx:`PredictionStrategy::DEFAULT_PREDICTION`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step.

.. warning::

   In mechanics, this may lead to very high increments of the deformation
   gradients or the strain in the neighboring elements of boundaries where
   evolving displacements are imposed.
   
Prediction for the state at the beginning of the time step (:cxx:`PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

Prediction based on a behaviour integration with constant gradients (:cxx:`PredictionStrategy::CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

FBar formulation
----------------

 The :math:`\bar{F}` method, or FBar formulation, is implemented following
 :cite:`de_souza_neto_design_1996`. This formulation is designed to handle
 nearly incompressible materials in large strain analysis by using a modified
 deformation gradient :math:`\bar{\underline{F}}` that separates volumetric and
 deviatoric responses.

 The method replaces the standard deformation gradient :math:`\underline{F}`
 with an assumed modified counterpart :math:`\bar{\underline{F}}` in the
 computation of stresses. This modification is based on a multiplicative split
 into volumetric and deviatoric parts:

 .. math::

    \underline{F} = J^{1/3}\, \bar{\underline{F}}

 where :math:`J = \det{\underline{F}}` is the Jacobian of the deformation
 gradient. This formulation effectively avoids locking issues in nearly
 incompressible materials while maintaining accuracy.

 The :math:`\bar{F}` formulation is particularly suited for low-order finite
 elements and is applicable to arbitrary material models. It ensures quadratic
 rates of convergence in Newton-Raphson schemes.

 This formulation is enabled by passing an additional parameter to the
 :cxx:`Mechanics` behaviour integrator, as follows:

 .. code:: c++

   const auto fbar_parameters = mfem_mgis::dict{
       {"Regularization", mfem_mgis::dict{{{"FBar", mfem_mgis::list{}}}}};
   mechanics.addBehaviourIntegrator(ctx, "Mechanics", library,
                                    behaviour, fbar_parameters) | or_die;

Faltus 2026 regularization
---------------------------

 The regularization proposed by Faltus et al. in the context of contact
 mechanics using a third medium is implemented here :cite:`faltus_deformation_2026`. This
 regularization only applies to finite strain behaviours. Currently, this
 regularization is only available for isotropic behaviours.

 This regularization adds a contribution to the standard variational
 operator in finite strain and can be derived from an energy :math:`W`
 which penalizes the difference between the deformation gradient
 :math:`\underline{F}` at a given quadrature point and its value
 :math:`\bar{\underline{F}}` at the centroid of the element:

 .. math::

     W\left(\underline{F}, \bar{\underline{F}}\right) =
     \alpha\,\left(\underline{F}-\bar{\underline{F}}\right)\,\colon\,
     \left(\underline{F}-\bar{\underline{F}}\right)

 where :math:`\alpha` is a penalization coefficient.

 This regularization is enabled by passing an additional parameter to the
 :cxx:`Mechanics` behaviour integrator, as follows:

 .. code:: c++

   const auto faltus_parameters = mfem_mgis::Parameters{
       {"Regularization",
        mfem_mgis::Parameters{
            {"Faltus2026",
             mfem_mgis::Parameters{{{"PenalizationCoefficient", 1e11}}}}}};
   mechanics.addBehaviourIntegrator(ctx, "Mechanics", "ThirdMedium", library,
                                    behaviour2, faltus_parameters) | or_die;

The :cxx:`info` function
------------------------

The :cxx:`info` function allows displaying information about an object
in an output stream.


Retrieving information on a finite element discretization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of usage
""""""""""""""""

.. code:: c++

   const auto& fed = problem.getFiniteElementDiscretization();
   const auto success = mfem_mgis::info(ctx, fed);

Retrieving information on a partial quadrature space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cxx:`PartialQuadratureSpaceInformation` structure contains some
relevant information about a partial quadrature space:

- the identifier and the name of the underlying material,
- the total number of elements,
- the total number of integration points,
- the number of elements points per geometric type.
- the number of quadrature points per geometric type.

This structure is created by:

- :cxx:`getLocalInformation`, which returns the information relative to
  the current process.
- :cxx:`getInformation`, which returns the information gathered from all
  processes.

The :cxx:`PartialQuadratureSpaceInformation` structure can be printed to
an output stream using the :cxx:`info` function.

Example of usage
""""""""""""""""

.. code:: c++

   const auto& qspace =
      problem.getBehaviourIntegrator(1).getPartialQuadratureSpace();
   const auto success = mfem_mgis::info(ctx, std::cout, qspace);

Evaluators of quantities at integration points
----------------------------------------------

Overview
^^^^^^^^

The :cxx:`setMaterialProperty` and :cxx:`setExternalStateVariable` methods of behaviour integrators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Resolution of dependencies of a nonlinear evolution problemn using another
--------------------------------------------------------------------------

The :cxx:`Simulation` class
---------------------------

Physical system, coupling schemes and models
--------------------------------------------

Line-search-like handling of behaviour integration failures
-----------------------------------------------------------

Profiling Toolkit
-----------------

Timers are now accessible by MGIS contexts, offering greater flexibility:

- Eliminates the need for a single global object
- Allows multiple contexts and the ability to distinguish different parts of the code managed through separate objects
- Prevents potential negative interactions with other code using the same timer type.

Example of usage:

.. code:: c++

    CatchTimeSection(ctx, "Class::FunctionName");

Refactoring of submesh and finite element space creation
--------------------------------------------------------

The creation of submeshes and associated finite element spaces has been refactored to be handled through the :cxx:`MeshDiscretization` and :cxx:`FiniteElementSpacesManager` classes. This refactoring introduces the following new classes:

- :cxx:`MeshDiscretization`: A fundamental class that handles the lifetime of meshes and provides utilities for managing mesh-related operations, including submesh creation based on material or boundary identifiers.
- :cxx:`FiniteElementDiscretization`: Extends :cxx:`MeshDiscretization` to handle the lifetime of finite element collections and spaces, providing a high-level interface for creating and managing finite element spaces.
- :cxx:`FiniteElementSpacesManager`: Manages similar finite element spaces (siblings) that share the same mesh and finite element collection but may have different vectorial dimensions, with automatic reuse of existing spaces.

This refactoring improves code organization, reduces duplication, and provides a more consistent and flexible API for working with submeshes and finite element spaces.

Search of the behaviour integrator providing a quantity
-------------------------------------------------------

The functions :cxx:`hasGradientProvider`,
:cxx:`hasThermodynamicForceProvider` and
:cxx:`hasInternalStateVariableProvider` search, among the behaviour
integrators defined on a given location, the one providing a given
quantity. The ``ParaviewExportIntegrationPointResultsAtNodes``
post-processing relies on them, so that a result can be exported on a
material carrying several behaviour integrators.

Issues fixed
============

- Issue 292: Remove deprecated usage of `getMaterial` in
  `ParaviewExportIntegrationPointResultsAtNodesBase::getPartialQuadratureFunctionViews`
  and `ParaviewExportIntegrationPointResultsAtNodesBase::getResultDescription`
- Issue 290: Add support for partial quadratures spaces on boundaries
- Issue 288: ￼master doesnt build without MPI
- Issue 286: updateGridFunction dilutes interface values when the
  functions do not cover all the materials of the mesh
- Issue 284: NonLinearModel does not export the initial state
- Issue 281:￼ Nodal export of integration point results uses the node
  index as an integration point index
- Issue 279: ParaviewExportResults silently ignores unknown parameters
- Issue 277: NonLinearModel calls the deprecated solve and loses the
  caller's context
- Issue 275: UniformHeatSourceBoundaryCondition: only the last
  integration point contributes to the residual ￼
- Issue 273: Fix forgotten change while parallel execution bug in
  behaviour integration ￼
- Issue 270: FirstIterationConvergenceCriterion never triggers coupling
  iterations + wrong post-processing time
- Issue 268: Parallel execution bug in behaviour integration
- Issue 264: Installation from the doc fails on a fresh setup
- Issue 262: Refactor creation of submeshes and associated finite
  element spaces to handle it through MeshDiscretization and
  FiniteElementSpacesManager enhancement ￼
- Issue 260: Refactor commented examples documentation
- Issue 257: Allow NewtonSolver to discard linear solver failures
- Issue 254: Missing parameter option for GMRESSolver
- Issue 253: Incomplete linear solver convergence checks in
  NonLinearEvolutionProblemImplementation.cxx and NewtonSolver.cxx
- Issue 248: Improve `PartialQuadratureFunction` interface
- Issue 245: Check the size of the unknown when defining bricks
- Issue 240: Small bug in `LinearSolverFactory.cxx`
- Issue 237: [cmake] Add a build-tests target
- Issue 218: [performance] synchronize success of the setup methods at a
  higher level to minimize collective communications enhancement
- Issue 213: ￼ Add a simple way to resolve dependencies (material
  properties, external state variables) of a :cxx:`NonLinearEvolutionProblem`
  using the gradients, thermodynamic forces and internal state variables
  of another one.
- Issue 211: Add coupling schemes and models
- Issue 209: Introduce the :cxx:`Simulation` class￼
- Issue 206: Line-search-like handling of behaviour integration failures
- Issue 200: automatically assign materials and boundaries's names from
  |mfem|'s attributes ￼
- Issue 198: Add the ability to define multiple behaviour integrators on
  the same material
- Issue 193: Add support for other types of search operators for
  computing a prediction of the solution at the end of the time step
  enhancement
- Issue 192: Take external forces into account when computing the
  prediction of the solution
- Issue 188: retrieve information about a quadrature space
- Issue 149: work on the prediction of the solution
