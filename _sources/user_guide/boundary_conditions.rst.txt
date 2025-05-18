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

The :code:`UniformImposedPressure` boundary condition computes the
residual associated with the following variational operator:

.. math::

   -\int_{\partial\Omega_{r}} p\,\vec{N}\,\cdot\,\,\vec{u}^{\star}\,\mathrm{d}S

where:

- :math:`\partial\Omega_{r}` denotes the boundary in the reference
  configuration on which the pressure is applied,
- :math:`p` is the imposed pressure,
- :math:`\vec{N}` is normal to the boundary in the reference configuration,
- :math:`\vec{u}^{\star}` is a virtual displacement.

The following snippet imposes an uniform pressure, which evolves linearly
with time, on the boundary labeled :code:`3`:

.. code:: c++
   
   problem.addBoundaryCondition(
      std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
          problem.getFiniteElementDiscretizationPointer(), 3,
          [](const mfem_mgis::real t) noexcept { return 150e6 * t; }));

Heat transfer boundary conditions
=================================

The :code:`UniformHeatSource` boundary condition computes the
residual associated with the following variational operator:

.. math::

   \int_{\Omega_{r}} q\,T^{\star}\,\mathrm{d}V

where:

- :math:`\Omega_{r}` denotes the domain of interest in the reference
  configuration,
- :math:`q` is the heat source,
- :math:`T^{\star}` is a virtual temperature.

The following snippet imposes an uniform heat source, which evolves
linearly with time, on the material named :code:`fuel`:

.. code:: c++
   
   problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformHeatSourceBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(), "fuel",
            [&q](const auto t) { return q * t; }));
