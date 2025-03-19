.. _mfem_mgis_boundary_conditions:

===================
Boundary conditions
===================

.. contents::
    :depth: 3
    :local:

General purpose boundary conditions
===================================

The :code:`UniformDirichletBoundaryCondition` boundary condition
----------------------------------------------------------------

The :code:`UniformDirichletBoundaryCondition` allows to impose
a prescribed value on a boundary.

The following snipet sets to zero the first component of the unknowns
on the boundary labeled :code:`5`:

.. code:: c++

   problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(), 5, 0));

Mechanical boundary conditions
==============================

The :code:`UniformImposedPressure` boundary condition
-----------------------------------------------------

The following snippet imposes a uniform pressure, which evolves linearly
with time, on the boundary labeled :code:`3`:

.. code:: c++
   
   problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3,
          [](const mfem_mgis::real t) noexcept { return 150e6 * t; }));

.. Heat transfer boundary conditions
.. =================================
