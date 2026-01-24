.. _mfem_mgis_partial_quadrature_functions:

============================
Partial quadrature functions
============================

Partial quadrature functions map quadrature points to values. A partial
quadrature function is defined on the quadrature points associated with
a unique material.

Projection of a grid function to a partial quadrature function
==============================================================

The :cxx:`update` function allows projecting a grid function to a
partial quadrature function.

Example of usage
----------------

.. code:: c++
   
   const auto success = update(ctx, dest, src);

Extrapolating quadrature functions at nodes
===========================================

Using a so-called :cxx:`GridFunction`, a set of partial quadrature
functions, each associated with a distinct material, can be extrapolated
at nodes:

- the :cxx:`makeGridFunction` creates the :cxx:`GridFunction` and the associated
  finite element space.
- the :cxx:`updateGridFunction` fills the :cxx:`GridFunction` with the values of
  the partial quadrature functions.

A :cxx:`GridFunction` being an :cxx:`MFEM` object, the reader may refer
to the :cxx:`MFEM` documentation for further details. In particular, a
:cxx:`GridFunction` can easily be exported to `Paraview` using
:cxx:`MFEM`'s :cxx:`ParaViewDataCollection`.

Projecting quadrature functions at nodes
========================================

The :math:`L_{2}` projection :math:`\bar{f}` of a scalar quadrature
function :math:`f` is obtained by minimizing the following functional:

.. math::

   \int_{\Omega} \left(\bar{f}-f\right)^{2}\,\mathrm{d}\,V

where :math:`\Omega` is the domain on which `f` is defined.

For a set of quadrature functions, the projection is computed on the
union of the domains of the quadrature functions.

For vectorial functions, the projection is computed component-wise.

Available functions
-------------------

The :cxx:`computeL2Projection` function allows projecting a set of
quadrature functions on nodes by minimizing the :math:`L_{2}` norm of
the difference of the projection and the values of the quadrature functions.

The :cxx:`computeL2Projection` function internally calls successively the
functions :cxx:`createL2ProjectionResult` and :cxx:`updateL2Projection`.
Those functions can be used directly if the projection has to be
performed several times.

Example of usage
----------------

.. code-block:: c++

   auto f = LinearFactor<parallel>::getFactory();
   auto fct = PartialQuadratureFunction::evaluate(
       fespace, [](const real x, const real y) noexcept { return cos(x) * y; });
   auto os = f.generate(ctx, "MUMPSSolver", fespace, {});
   const auto oresult = computeL2Projection<parallel>(ctx, *os, {*fct});

Implicit gradient regularization
================================

For a scalar function :math:`f`, the implicit gradient regularization
:math:`\bar{f}` is defined as the solution of:

.. math::

   \bar{f}-l_{c}^{2}\cdot\Delta\bar{f}=f

where :math:`l_{c}` is a characteristic length.

See :cite:`peerlings_gradient_1996` for details.

For a vectorial function, the implicit gradient regularization is
computed component-wise.

Available functions
-------------------

The :cxx:`computeImplicitGradientRegularization` function solves
the previous equation to determine :math:`\bar{f}`.

The :cxx:`computeImplicitGradientRegularization` function internally calls
successively the functions :cxx:`createL2Projection` and
:cxx:`updateImplicitGradientRegularization`. Those functions can be used
directly if the projection has to be performed several times.

Example of usage
----------------

.. code-block:: c++

   auto fct = PartialQuadratureFunction::evaluate(
       space, [](const real x, const real) noexcept {
         const auto c = std::cos(4 * pi * x);
         return (1 + l0 * l0 * 16 * pi * pi) * c;
       });
   auto os = f.generate(ctx, "MUMPSSolver", fespace, {});
   const auto oresult =
       computeImplicitGradientRegularization<parallel>(ctx, *os, {*fct}, l0);


Miscellaneous functions
=======================

:cxx:`rotateThermodynamicForcesToGlobalFrame`
---------------------------------------------

