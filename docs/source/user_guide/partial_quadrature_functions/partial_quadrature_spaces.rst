.. _mfem_mgis_partial_quadrature_spaces:

=========================
Partial quadrature spaces
=========================

Retrieving information on a partial quadrature space
====================================================

The :cxx:`PartialQuadratureSpaceInformation` structure contains some
relevant information about a partial quadrature space:

- the identifier and the name of the underlying material,
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

   const auto& qspace =
      problem.getBehaviourIntegrator(1).getPartialQuadratureSpace();
   const auto success = mfem_mgis::info(ctx, std::cout, qspace);
