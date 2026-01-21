.. _mfem_mgis_partial_quadrature_functions:

============================
Partial quadrature functions
============================

Partial quadrature functions map quadrature points to values. A partial
quadrature function is defined on the quadrature points associated with
a unique material.

Extrapolating quadrature functions at nodes
===========================================

Using a so-called :cxx:`GridFunction`, a set of partial quadrature
functions, each associated with a distinct material, can be extrapolated
at nodes:

- the `makeGridFunction` creates the `GridFunction` and the associated
  finite element space.
- the `updateGridFunction` fills the `GridFunction` with the values of
  the partial quadrature functions.

A :cxx:`GridFunction` being an `MFEM` object, the reader may refer to
the `MFEM` documentation for further details. In particular, a
:cxx:`GridFunction` can easily be exported to `Paraview` using
:cxx:`MFEM`'s :cxx:`ParaViewDataCollection`.

Projecting quadrature functions at nodes
========================================

The :cxx:`computeL2Projection` function allows projecting a set of
quadrature functions on nodes by minimizing the L2 norm of the
projection and the values of the quadrature functions.

The :cxx:`computeL2Projection` internally calls successively the
functions :cxx:`createL2ProjectionResult` and `updateL2Projection`.
Those functions can be used directly if the projection has to be
performed several times.

Available functions
===================

:cxx:`rotateThermodynamicForcesToGlobalFrame`
---------------------------------------------

