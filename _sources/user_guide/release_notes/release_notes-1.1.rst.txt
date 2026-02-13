.. _mfem_mgis_release_notes_1_1:

============================
Release notes of Version 1.1
============================

Documentation
=============

Doxygen documentation
---------------------

The `doxygen` documentation is available on `this page <https://thelfer.github.io/mfem-mgis-doxygen/index.html>`_.

New features
============

Projection of a grid function to a partial quadrature function
--------------------------------------------------------------

The :cxx:`update` function allows projecting a grid function to a
partial quadrature function.


:math:`L_{2}` projection of a partial quadrature function to a grid function
----------------------------------------------------------------------------

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
  :cxx:`NonLinearEvolutionProblemImplementation<false>`. See Issue 168
  for details.
- The :cxx:`update` allows to update the values for a
  `PartialQuadratureFunction` from the values of a `GridFunction`. See
  Issue 169 for details.

Deprecated functions and methods
--------------------------------

- :cxx:`getUnknownsAtTheBeginningOfTheTimeStep` and
  :cxx:`getUnknownsAtTheEndOfTheTimeStep` are deprecated in favor of
  :cxx:`getUnknowns` in :cxx:`AbstractNonLinearEvolutionProblem`,
  :cxx:`NonLinearEvolutionProblem`, and
  :cxx:`NonLinearEvolutionProblemImplementationBase`.

Issues fixed
============

- Issue 178: Deployment of the `doxygen` documentation
- Issue 176: Modify `MFEMMGISConfig.cmake.in` to define `TFEL_DIR` and
  `MFrontGenericInterface_DIR` using by default the values passed to
  `cmake` at the configuration stage
- Issue 175: Update `tests/README.md`
- Issue 169: Add the ability to project a :cxx:`GridFunction` to a
  :cxx:`PartialQuadratureFunction`.
- Issue 168: Add helper methods to wrap the unknowns of a non linear
  evolution problem as :cxx:`GridFunction`.
- Issue 167: Allow tensorial external state variables
- Issue 162: Change the API of :cxx:`LinearSolverFactory`
- Issue 159: [fix] :cxx:`upgradeGridFunction` seems wrong for grid
  functions defined on a :cxx:`SubMesh` ￼ Status: Closed (completed).
- Issue 152: [post-processing] Add support for projecting ip fields at
  nodes using L2 projection
