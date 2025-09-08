===========
Mesh Reader
===========

``Mfem-Mgis`` uses the classic sequential reading of meshes generated via gmsh by specifying the `MFEM` input format with the ``-format msh22`` option. 

.. note::
  For more information: https://mfem.org/mesh-formats/

Keywords
^^^^^^^^

+------------------+-------------------------------------------------------------------+
| Key              | Description                                                       |
+==================+===================================================================+
| ``MeshFileName`` | Define the mesh file name.                                        |
+------------------+-------------------------------------------------------------------+
| ``MeshReadMode`` | Read Splitted Mesh, FromScratch = sequential, Restart = parallel  |
+------------------+-------------------------------------------------------------------+
 
.. code-block::

   mfem_mgis::Parameters{{"MeshFileName", "mesh-explorer.mesh."},
                         {"MeshReadMode", "Restart"}, ... }

Split Mesh
^^^^^^^^^^

If your mesh is too large (memory limit), you can also use MFEM's option of spliting the mesh into smaller meshes and reading them in parallel.

To split the mesh, we use the ``mesh-explorer`` tool (https://mfem.org/meshing-miniapps/#mesh-explorer) and you can specify the parallel reader using the keywork: ``MeshReadMode``.


Install Mesh-explorer
---------------------

Here are two ways to install mesh-explorer. We recommend installing via spack, but we will also add the installation procedure with cmake.

Spack version:

.. code-block::

  spack install mfem+miniapps  

The binary file is then located at: 

.. code-block::

  `spack location -i mfem`/share/mfem/miniapps/meshing/mesh-explorer

.. note::

  If you installed mfem-mgis using the spack installation procedure. The installation has already been performed with the miniapps.


CMake version:

.. code-block::

   git clone https://github.com/mfem/mfem.git
   mkdir build && cd build
   spack load hypre
   cmake ../mfem/ -DMFEM_USE_MPI=ON -DMFEM_ENABLE_MINIAPPS=ON
   make -j 8

The binary file is then located at: 

.. code-block::

  build/miniapps/meshing/mesh-explorer


Mesh generation
---------------

The mesh explorer works interactively. Here is a series of options for splitting a gmsh mesh named into `128` meshes named `output-mesh.XXXXXX`

Run:

.. code-block::

  mesh-explorer --mesh mesh.msh


Create partitioning:

.. code-block::

  press p > for Generate a partitioning
  press 1 > Choose the partitioning method
  press 128 > choose the number of MPI processes
  
Save mesh in parallel format:

.. code-block::

  press T > choose to save in parallel format
  press output-mesh. > prefix name for output files
  press 6 > number of digit


Load your mesh with MFEM-MGIS
-----------------------------

.. note::

  Please, verify that the splited mesh chosen is splited in 128 files if you want to run it on 128 mpi processes.

Example: 

- One example is available here: tests/ParallelReadMode.cxx

.. code-block:: c++

    auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{{"MeshFileName", "output-mesh."},
                            {"MeshReadMode", "Restart"},
                            {"FiniteElementFamily", "H1"},
                            {"FiniteElementOrder", p.order},
                            {"UnknownsSize", dim},
                            {"NumberOfUniformRefinements", p.parallel ? 2 : 0},
                            {"Parallel", p.parallel}});
