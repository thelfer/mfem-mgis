============
Getting help
============

The :cxx:`info` function
========================

The :cxx:`info` function allows to display information about an object
in an output stream.

Example of usage
----------------

.. code:: c++

   const auto& fed = problem.getFiniteElementDiscretization();
   const auto success = mfem_mgis::info(ctx, fed);
