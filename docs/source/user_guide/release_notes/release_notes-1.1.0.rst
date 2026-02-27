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

- The name of the materials and boundaries are automatically retrieved
  from |MFEM|'s mesh.
- Many methods have been deprecated to have a consistent error handling
  scheme based on `MGIS`'s one. As such, many methods and functions now
  takes and `MGIS`'s :cxx:`Context` as their first argument.

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
- :math:`\Delta\,\mathbb{u}` is submitted to the the increment of the
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
  gradients are constant over the time (and thus egal to their values at
  the beginning of the time step).
- :math:`\Delta\,\mathbb{u}` is submitted to the the increment of the
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
- :cxx:`IntegrationOperator::CONSISTENT_TANGENT`: the constistent
  tangent operator, defined by the derivative of the thermodynamic force
  with respect to the gradients at the end of the time step. See
  :cite:`simo_consistent_1985` for details.

The :cxx:`info` function
------------------------

The :cxx:`info` function allows to display information about an object
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

Issues fixed
============

- Issue 200: automatically assign materials and boundaries's names from
  MFEM's attributes ￼
- Issue 188: retrieve informations about a quadrature space
- Issue 149: work on the prediction of the solution
