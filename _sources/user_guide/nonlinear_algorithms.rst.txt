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

   mechanics.setPredictionPolicy(
      {.strategy = PredictionStrategy::ELASTIC_PREDICTION});

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
   
Elastic prediction (:cxx:`PredictionStrategy::ELASTIC_PREDICTION`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :cxx:`ElasticPrediction` strategy determines the increment of the
displacement :math:`\Delta\,\mathbb{u}` by solving the following linear
system:

.. math::

   \mathbb{K}_{e}\,\cdot\,\Delta\,\mathbb{u} = \ets{\mathbb{F}_{e}}-\bts{\mathbb{F}_{i}}

where:

- :math:`\mathbb{K}_{e}` denotes the elastic stiffness matrix.
- :math:`\bts{\mathbb{F}_{e}}` denotes the external forces at the beginning
  of the time step.
- :math:`\bts{\mathbb{F}_{i}}` denotes the inner forces at the beginning
  of the time step.
- :math:`\Delta\,\mathbb{u}` is submitted to the the increment of the
  imposed Dirichlet boundary conditions.

.. note::

   Although the wording explicitly refers to mechanics, this equation
   applies to all physics.
