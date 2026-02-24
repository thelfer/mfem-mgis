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

   mechanics.setPredictionPolicy(
      {.strategy = PredictionStrategy::ELASTIC_PREDICTION});

Available strategies
^^^^^^^^^^^^^^^^^^^^

Default prediction
""""""""""""""""""

By default, a nonlinear evolution problem uses the solution at the
beginning of the time step, modified by applying Dirichlet boundary
conditions, as the initial guess of the solution at the end of the time
step.

.. warning::

   In mechanics, this may lead to very high increments of the deformation
   gradients or the strain in the neighboring elements of boundaries where
   evolving displacements are imposed.
   
Elastic prediction
""""""""""""""""""

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

Retrieving information on a partial quadrature space
----------------------------------------------------

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
^^^^^^^^^^^^^^^^

.. code:: c++

   const auto& qspace =
      problem.getBehaviourIntegrator(1).getPartialQuadratureSpace();
   const auto success = mfem_mgis::info(ctx, std::cout, qspace);

Issues fixed
============

- Issue 
- Issue 149: work on the prediction of the solution