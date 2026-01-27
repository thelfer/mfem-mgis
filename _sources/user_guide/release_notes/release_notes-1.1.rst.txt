.. _mfem_mgis_release_notes_1_1:

============================
Release notes of Version 1.1
============================

New features
============

The :math:`L_{2}` projection :math:`\bar{f}` of a scalar quadrature
function :math:`f` is obtained by minimizing the following functional:

.. math::

   \int_{\Omega} \left(\bar{f}-f\right)^{2}\,\mathrm{d}\,V

where :math:`\Omega` is the domain on which `f` is defined.

For a set of quadrature functions, the projection is computed on the
union of the domains of the quadrature functions.

For vectorial functions, the projection is computed component-wise.

The :cxx:`computeL2Projection` function allows projecting a set of
quadrature functions on nodes by minimizing the :math:`L_{2}` norm of
the difference of the projection and the values of the quadrature
functions.

The :cxx:`computeL2Projection` internally function calls successively
the functions :cxx:`createL2ProjectionResult` and
:cxx:`updateL2Projection`. Those functions can be used directly if the
projection has to be performed several times.

Implicit gradient regularization
--------------------------------

For a scalar function :math:`f`, the implicit gradient regularization
:math:`\bar{f}` is defined as the solution of:

.. math::

   \bar{f}-l_{c}^{2}\cdot\Delta\bar{f}=f

where :math:`l_{c}` is a characteristic length.

See :cite:`peerlings_gradient_1996` for details.

For a vectorial function, the implicit gradient regularization is
computed component-wise.

The :cxx:`computeImplicitGradientRegularization` function allows
projecting a set of quadrature functions on nodes by solving the
previous equation.

The :cxx:`computeImplicitGradientRegularization` function internally
calls successively the functions :cxx:`createL2Projection` and
:cxx:`updateImplicitGradientRegularization`. Those functions can be used
directly if the projection has to be performed several times.

Application programming interface
=================================

New functions and methods
-------------------------

- :cxx:`getUnknownsAsGridFunction` has been added to
  :cxx:`NonLinearEvolutionProblemImplementation<true>` and
  :cxx:`NonLinearEvolutionProblemImplementation<false>`.

Deprecated functions and methods
--------------------------------

- :cxx:`getUnknownsAtTheBeginningOfTheTimeStep` and
  :cxx:`getUnknownsAtTheEndOfTheTimeStep` are deprecated in favor of
  :cxx:`getUnknowns` in :cxx:`AbstractNonLinearEvolutionProblem`,
  :cxx:`NonLinearEvolutionProblem`, and
  :cxx:`NonLinearEvolutionProblemImplementationBase`.

Issues fixed
============

- Issue 167: Add helper methods to wrap the unknowns of a non linear
  evolution problem as :cxx:`GridFunctions`
- Issue 167: Allow tensorial external state variables
- Issue 162: Change the API of :cxx:`LinearSolverFactory`
- Issue 159: [fix] :cxx:`upgradeGridFunction` seems wrong for grid
  functions defined on a :cxx:`SubMesh` ￼ Status: Closed (completed).
- Issue 152: [post-processing] Add support for projecting ip fields at
  nodes using L2 projection
