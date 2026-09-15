.. _mfem_mgis_MeshSetDiscretization:

=======================
MeshSetDiscretization
=======================

The ``MeshDiscretization`` class is a fundamental component of ``MFEM-MGIS`` that handles the lifetime of meshes and provides utilities for managing mesh-related operations. It serves as the base class for more specialized discretization classes.

Main Features
-------------

- **Mesh Management**: Handles the creation, storage, and access to both sequential and parallel meshes.
- **Submesh Creation**: Provides methods to create submeshes based on material or boundary identifiers.
- **Identifier Management**: Supports mapping between material/boundary names and their numeric identifiers using regular expressions.
- **Parallel and Sequential Support**: Works with both parallel (MPI-based) and sequential meshes, determined by the ``Parallel`` parameter.
- **Point and Point Set Management**: Allows registration and retrieval of points and point sets in 2D and 3D.

Key Methods
-----------

- ``getMesh<parallel>()``: Returns the underlying mesh (parallel or sequential).
- ``getSubMesh<parallel>()``: Creates or retrieves a submesh based on material or boundary identifiers.
- ``getMaterialsIdentifiers()`` / ``getBoundariesIdentifiers()``: Resolves identifiers from names or numeric values.
- ``isDefinedOnMaterials()`` / ``isDefinedOnBoundaries()``: Checks if a mesh is defined on specific materials or boundaries.
- ``manages()``: Checks if a given mesh is managed by this discretization.

Usage Example
-------------

The following example demonstrates basic usage of ``MeshDiscretization``:

.. code-block:: c++

  using namespace mfem_mgis;
  auto ctx = mgis::Context{};
  auto omesh = construct<MeshDiscretization>(
      ctx,
      Parameters{{{"MeshFileName", "mesh.mesh"}},
                 {"NumberOfUniformRefinements", 0},
                 {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                 {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                 {"Parallel", false}});
  
  // Get a submesh for material "Attr1"
  const auto osm1 = omesh->getSubMesh<false>(
      ctx, "Attr1", MeshDiscretization::Location::ON_MATERIALS);
  
  // Check if a mesh is defined on materials
  const auto ook = omesh->isDefinedOnMaterials(ctx, *osm1);

Parameters
----------

- ``MeshFileName`` (string): Path to the mesh file.
- ``NumberOfUniformRefinements`` (int): Number of uniform refinements to apply to the mesh.
- ``Materials`` (dictionary): Mapping between material names and their identifiers.
- ``Boundaries`` (dictionary): Mapping between boundary names and their identifiers.
- ``Parallel`` (boolean): Whether to use parallel (MPI-based) computation.
- ``GeneralVerbosityLevel`` (int): Verbosity level for debugging output.

.. note::

  The ``MeshDiscretization`` class is designed to be lightweight, movable, and copyable. It uses the PIMPL idiom for its internal implementation.
