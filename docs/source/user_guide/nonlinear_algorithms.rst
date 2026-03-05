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

This reduction of the amplitude repeated until a estimation of the
solution that does not lead to an integration failure