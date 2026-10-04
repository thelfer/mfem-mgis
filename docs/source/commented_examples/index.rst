==================
Commented Examples
==================

The examples come from two repositories.

mfem-mgis-examples
==================

These examples are in https://github.com/latug0/mfem-mgis-examples.

``ctest`` runs their tests. The CMake option ``MFEM_MGIS_EXAMPLES_TEST_MODE``
selects their size:

- ``full`` runs the complete simulations.
- ``restricted`` only runs their beginning or a smaller case. It is much
  faster.
- ``auto``, the default, selects ``restricted`` in the ``Debug`` and
  ``Coverage`` builds and ``full`` otherwise.

Some tests are the same in both modes.

.. toctree::
   :maxdepth: 1

   basic_tests
   rve_elastic_inclusions
   rve_mox
   rjh_plate

mm-opera-hpc
============

These examples are in https://github.com/rprat-pro/mm-opera-hpc. They were
developed in the OperaHPC project.

.. toctree::
   :maxdepth: 1

   bubbles
   polycrystal
   cermet
