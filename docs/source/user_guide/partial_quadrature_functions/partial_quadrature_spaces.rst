.. _mfem_mgis_partial_quadrature_spaces:

=========================
Partial quadrature spaces
=========================

Creating a partial quadrature space
===================================

A partial quadrature space is defined by a finite element
discretization, a location, i.e. a material or a boundary, and a
function selecting the integration rule of each element:

.. code:: c++

   auto qspace = std::make_shared<PartialQuadratureSpace>(
       fed, BoundaryIdentifier{.id = 2}, selector);

A partial quadrature space defined on a boundary is built on the
boundary elements of the mesh of the finite element discretization.

Another constructor takes a finite element space as second argument, for
instance a finite element space defined on a submesh. This finite
element space must be managed by the finite element discretization.

Retrieving information on a partial quadrature space
====================================================

The :cxx:`PartialQuadratureSpaceInformation` structure contains some
relevant information about a partial quadrature space:

- the identifier and a description of the underlying material or
  boundary, including its name, if defined,
- the total number of elements,
- the total number of integration points,
- the number of quadrature points per geometric type.

This structure is created by:

- :cxx:`getLocalInformation`, which returns the information relative to
  the current process.
- :cxx:`getInformation`, which returns the information gathered from all
  processes.

The :cxx:`PartialQuadratureSpaceInformation` structure can be printed to
an output stream using the :cxx:`info` function.

.. note::

   There is also an overload of the :cxx:`info` function that
   takes a partial quadrature space as argument, which builds
   internally a :cxx:`PartialQuadratureSpaceInformation`,
   synchronizes it over all processes and prints the result.

Example of usage
----------------

.. code:: c++

   auto& bi = problem.getBehaviourIntegrator(ctx, 1, 0) | or_die;
   const auto& qspace = bi.getPartialQuadratureSpace();
   const auto success = mfem_mgis::info(ctx, std::cout, qspace);
