.. _mfem_mgis_FiniteElementSpaces:

========================
FiniteElementSpaces
========================

The ``FiniteElementSpacesManager`` class manages similar finite element spaces (siblings) that share the same mesh and finite element collection but may have different vectorial dimensions. It provides a centralized way to create, reuse, and manage finite element spaces.

Main Features
-------------

- **Finite Element Space Creation**: Creates finite element spaces on the whole mesh or on submeshes.
- **Space Reusability**: Reuses existing finite element spaces when the same configuration is requested.
- **Submesh Support**: Creates finite element spaces on submeshes defined by material or boundary identifiers.
- **Lightweight Design**: Designed to be lightweight, movable, and copyable.
- **Integration with MeshDiscretization**: Works seamlessly with ``MeshDiscretization`` to manage mesh-related operations.

Key Methods
-----------

- ``getFiniteElementSpace<parallel>(ctx, nc)``: Creates or reuses a finite element space with a given vectorial dimension.
- ``getFiniteElementSpace<parallel>(ctx, args)``: Creates or reuses a finite element space on a submesh using ``GetFiniteElementSpaceOnSubMeshArguments``.
- ``getFiniteElementCollection()``: Returns the underlying finite element collection.
- ``getMeshDiscretization()``: Returns the associated mesh discretization.
- ``manages()``: Checks if a given finite element space is managed by this manager.
- ``setNodalFiniteElementSpace()``: Assigns a suitable nodal finite element space to the underlying mesh.

Usage Example
-------------

The following example demonstrates basic usage of ``FiniteElementSpacesManager``:

.. code-block:: c++

  using namespace mfem_mgis;
  auto ctx = mgis::Context{};
  
  // Create a FiniteElementSpacesManager from parameters
  auto om = construct<FiniteElementSpacesManager>(
      ctx, dict{{"MeshFileName", "mesh.mesh"},
                {"FiniteElementFamily", "H1"},
                {"FiniteElementOrder", 2},
                {"NumberOfUniformRefinements", 0},
                {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                {"Parallel", false}});
  
  // Create a finite element space on material "Attr1"
  auto fes1 = om->getFiniteElementSpace<false>(
      ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
            .identifiers = list{"Attr1"},
            .number_of_components = 3});
  
  // The same space will be reused if requested again
  auto fes2 = om->getFiniteElementSpace<false>(
      ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
            .identifiers = 1,
            .number_of_components = 3});
  
  // fes1 and fes2 point to the same finite element space
  assert(fes1.get() == fes2.get());

Parameters
----------

- ``MeshFileName`` (string): Path to the mesh file.
- ``FiniteElementFamily`` (string): Name of the finite element family (e.g., ``H1``).
- ``FiniteElementOrder`` (int): Order of the polynomial approximation.
- ``NumberOfUniformRefinements`` (int): Number of uniform refinements to apply to the mesh.
- ``Materials`` (dictionary): Mapping between material names and their identifiers.
- ``Boundaries`` (dictionary): Mapping between boundary names and their identifiers.
- ``Parallel`` (boolean): Whether to use parallel (MPI-based) computation.

.. note::

  The ``FiniteElementSpacesManager`` class is particularly useful when working with multiple finite element spaces that share the same mesh and finite element collection, as it automatically handles reuse and memory management.
