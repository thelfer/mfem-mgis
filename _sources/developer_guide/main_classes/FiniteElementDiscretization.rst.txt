.. _mfem_mgis_FiniteElementDiscretization:

===============================
FiniteElementDiscretization
===============================

The ``FiniteElementDiscretization`` class extends ``MeshDiscretization`` to handle the lifetime of finite element collections and spaces. It provides a high-level interface for creating and managing finite element spaces on meshes or submeshes.

Main Features
-------------

- **Finite Element Space Management**: Creates and manages finite element spaces for both sequential and parallel computations.
- **Finite Element Collection Handling**: Stores and provides access to finite element collections.
- **Submesh Support**: Allows creation of finite element spaces on submeshes defined by material or boundary identifiers.
- **Reusability**: Reuses existing finite element spaces when possible to optimize memory and performance.
- **Integration with MeshDiscretization**: Inherits all mesh management capabilities from ``MeshDiscretization``.

Key Methods
-----------

- ``getFiniteElementSpace<parallel>()``: Returns the finite element space (parallel or sequential).
- ``getFiniteElementCollection()``: Returns the underlying finite element collection.
- ``getFiniteElementSpacesManager()``: Returns the manager for finite element spaces.
- ``getVSize()`` / ``getTrueVSize()``: Returns the total number of unknowns.
- ``setNodalFiniteElementSpace()``: Assigns a suitable nodal finite element space to the underlying mesh.

Usage Example
-------------

The following example demonstrates basic usage of ``FiniteElementDiscretization``:

.. code-block:: c++

  using namespace mfem_mgis;
  auto ctx = mgis::Context{};
  auto ofed = construct<FiniteElementDiscretization>(
      ctx, dict{{"MeshFileName", "mesh.mesh"},
                {"FiniteElementFamily", "H1"},
                {"FiniteElementOrder", 2},
                {"UnknownsSize", 2},
                {"NumberOfUniformRefinements", 0},
                {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                {"Parallel", false}});
  
  // Get the finite element space
  auto& fes = ofed->getFiniteElementSpace<false>();
  
  // Get the finite element spaces manager
  auto fespaces = ofed->getFiniteElementSpacesManager();
  
  // Create a finite element space on a submesh
  auto fes1 = fespaces.getFiniteElementSpace<false>(
      ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
            .identifiers = list{"Attr1"},
            .number_of_components = 3});

Parameters
----------

In addition to the parameters supported by ``MeshDiscretization``, ``FiniteElementDiscretization`` supports:

- ``FiniteElementFamily`` (string): Name of the finite element family (e.g., ``H1``).
- ``FiniteElementOrder`` (int): Order of the polynomial approximation.
- ``UnknownsSize`` (int): Number of components of the unknowns.

.. note::

  The ``FiniteElementDiscretization`` class is designed to be lightweight and integrates seamlessly with the ``FiniteElementSpacesManager`` for advanced finite element space management.
