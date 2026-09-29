Basic Tests
===========

.. contents::

TensileTest
-----------

website : https://github.com/latug0/mfem-mgis-examples/tree/master/ex1

Description:

This example is a cyclic tension-compression test on a unit cube. The
imposed axial strain goes up to 0.9 %, down to -2.1 % and back up to
1.9 %, in 100 time steps. Its test compares the results to the reference
file ``Plasticity.ref``.

.. figure:: img/ex1Start.png
    :alt: Illustration of the start of the TensileTest simulation.

.. figure:: img/ex1End.png
    :alt: Illustration of the end of the TensileTest simulation.

Problem Solved
~~~~~~~~~~~~~~

.. code:: text

   Plastic behaviour with linear isotropic hardening (Plasticity.mfront):
   [ parameters      , material ]
   [ Young Modulus   , 70e9     ];
   [ Poisson Ratio   , 0.34     ];
   [ Yield Stress    , 300e6    ];
   [ Hardening Slope , 10e9     ];

   Boundary conditions:
   - symmetry on the faces x = 0, y = 0 and z = 0
   - imposed displacement along x on the face x = 1

   Solver : Conjugate Gradient (default)

   Element:
   - Family H1
   - Order 1

Run This Simulation
~~~~~~~~~~~~~~~~~~~

Without any option, the example runs as described above, in serial or
in parallel. The ``-r`` option compares the results to the reference file:

.. code-block:: bash

   ./UniaxialTensileTestEx
   mpirun -n 2 ./UniaxialTensileTestEx -p 1
   ./UniaxialTensileTestEx -r Plasticity.ref

Available options
~~~~~~~~~~~~~~~~~

To customize the simulation, several options are available, as detailed
below.

+---------------------------------+--------------------------------------------+
| Command line                    | Description                                |
+=================================+============================================+
| --mesh or -m                    | Mesh file (default = cube.mesh)            |
+---------------------------------+--------------------------------------------+
| --reference-file or -r          | Reference file, compared to the results    |
|                                 | when given (default = none)                |
+---------------------------------+--------------------------------------------+
| --behaviour or -b               | Name of the behaviour                      |
|                                 | (default = Plasticity)                     |
+---------------------------------+--------------------------------------------+
| --internal-state-variable or -v | Internal state variable compared to the    |
|                                 | reference                                  |
|                                 | (default = EquivalentPlasticStrain)        |
+---------------------------------+--------------------------------------------+
| --library or -l                 | Material library                           |
|                                 | (default = src/libBehaviour.so)            |
+---------------------------------+--------------------------------------------+
| --linearsolver or -ls           | Linear solver. Serial: 0 -> CG,            |
|                                 | 1 -> GMRES, 2 -> UMFPack. Parallel:        |
|                                 | 0 -> CG, 1 -> GMRES, 2 -> HypreFGMRES,     |
|                                 | 3 -> MUMPS (HyprePCG without MUMPS),       |
|                                 | 4 -> HypreGMRES (default = 0)              |
+---------------------------------+--------------------------------------------+
| --order or -o                   | Finite element order (polynomial degree)   |
|                                 | (default = 1)                              |
+---------------------------------+--------------------------------------------+
| --parallel or -p                | 0 for a serial run, 1 for a parallel run   |
|                                 | (default = 0)                              |
+---------------------------------+--------------------------------------------+

TwoLayerCube
------------

website: https://github.com/latug0/mfem-mgis-examples/tree/master/ex3

Description:

Periodic unit cube made of two elastic layers. A macroscopic strain is
imposed. The solution is compared to the analytical one.

Problem solved
~~~~~~~~~~~~~~

.. code:: text

   The layers are split at x = 0.5. Material 1 fills x < 0.5. Material 2
   fills x > 0.5.

   One component of the macroscopic strain is imposed:
   Exx -> 0, Eyy -> 1, Ezz -> 2, Exy -> 3, Exz -> 4, Eyz -> 5.
   Its value is 1 for a normal component. It is √2/2 for a shear
   component, in Mandel notation.

   Solver : CGSolver

   Elastic behaviour parameters, IsotropicLinearElasticity.mfront :
   [ parameters             , material 1 , material 2 ]
   [ First Lame Coefficient , 100        , 200        ];
   [ Shear Modulus          , 75         , 150        ];

   Element:
   - Family H1
   - Order 1

Run the simulation
~~~~~~~~~~~~~~~~~~

The first command runs the default case Eyy. The second one runs the case
Exy on 4 processes:

.. code-block:: bash

   ./two_layer_cube
   mpirun -n 4 ./two_layer_cube --test-case 3

The mesh ``cube_2mat_per.mesh`` has 4x4x4 hexahedra. The mesh ``Box.med``
has 8x8x8 hexahedra. Its periodicity is described by ``Box.per``. Reading it
requires MFEM built with MED support:

.. code-block:: bash

   ./two_layer_cube --mesh Box.med

Available options
~~~~~~~~~~~~~~~~~

+----------------------------------+----------------------------------+----------------------+
| Command line                     | Description                      | Default              |
+==================================+==================================+======================+
| ``--mesh`` or ``-m``             | Mesh file                        | cube_2mat_per.mesh   |
+----------------------------------+----------------------------------+----------------------+
| ``--library`` or ``-l``          | Material library                 | src/libBehaviour.so  |
+----------------------------------+----------------------------------+----------------------+
| ``--order`` or ``-o``            | Finite element order             | 1                    |
+----------------------------------+----------------------------------+----------------------+
| ``--xmax`` or ``-xm``,           | Coordinates of the upper corner  | 1                    |
| ``--ymax`` or ``-ym``,           | of the cube. They must match the |                      |
| ``--zmax`` or ``-zm``            | mesh.                            |                      |
+----------------------------------+----------------------------------+----------------------+
| ``--test-case`` or ``-t``        | Imposed component of the strain, | 1                    |
|                                  | from 0 to 5                      |                      |
+----------------------------------+----------------------------------+----------------------+
| ``--linearsolver`` or ``-ls``    | Linear solver: GMRESSolver,      | CGSolver             |
|                                  | CGSolver, UMFPackSolver or       |                      |
|                                  | MUMPSSolver. UMFPackSolver is    |                      |
|                                  | sequential only. MUMPSSolver is  |                      |
|                                  | parallel only.                   |                      |
+----------------------------------+----------------------------------+----------------------+
| ``--parallel`` or ``-p``,        | Run in parallel or not           | parallel             |
| ``--no-parallel`` or ``-no-p``   |                                  |                      |
+----------------------------------+----------------------------------+----------------------+
| ``--check`` or ``-c``,           | Compare or not the solution to   | compare              |
| ``--no-check`` or ``-no-c``      | the analytical one. It is only   |                      |
|                                  | valid for the provided meshes.   |                      |
+----------------------------------+----------------------------------+----------------------+

Satoh
-----

website: https://github.com/latug0/mfem-mgis-examples/tree/master/ex5

Description:

Plate of length 1 in plane strain, clamped on its left and right boundaries.
A parabolic temperature profile is imposed along the x-axis.

.. figure:: img/SatohTest.png
    :alt: Illustration of the displacement of the plate.

Problem solved
~~~~~~~~~~~~~~

.. code:: text

   This test models a 2D plate of length 1 in plane strain clamped on the left
   and right boundaries and subjected to a parabolic temperature profile along
   the x-axis:

   - the temperature is 293.15 K on the left and right boundaries
   - the temperature is 2000 K for x = 0.5

   This example shows how to define an external state variable using an
   analytical profile.

   Solver : UMFPackSolver
   Preconditioner : None

   Thermoelastic behavior law parameters :
   [ parameters            , material ]
   [ Young Modulus         , 150e9    ];
   [ Poisson Ratio         , 0.3      ];
   [ Thermal Expansion     , 1e-5     ];
   [ Reference Temperature , 293.15   ];

   Element: 
   - Family H1
   - Order 2

Run the simulation
~~~~~~~~~~~~~~~~~~

Parameters are hardcoded in this example.

.. code-block:: bash

   ./satoh

.. note::

   The example runs sequentially. A parallel run needs ``parallel`` set to
   ``true`` in the source code and a parallel linear solver.

Ssna303 Example (2D and 3D)
---------------------------

- website 2D example: https://github.com/latug0/mfem-mgis-examples/tree/master/ex2
- website 3D example: https://github.com/latug0/mfem-mgis-examples/tree/master/ex4

Description:

Tensile test on a notched beam with a finite-strain plastic behaviour. The 2D
example ex2 is in plane strain. The :ref:`tutorial <mfem_mgis_tutorial>`
describes it. The 3D example ex4 is described below.

.. figure:: img/ssna303Start.png
    :alt: Illustration of the start of the 3D ssna303 simulation.

.. figure:: img/ssna303End.png
    :alt: Illustration of the end of the 3D ssna303 simulation.

Problem solved in 3D
~~~~~~~~~~~~~~~~~~~~

.. code:: text

   The mesh ssna303_3d.msh is made of hexahedra. Lengths are in meters.
   The mesh is 1.5e-3 thick.

   Plastic behaviour with linear isotropic hardening,
   IsotropicLinearHardeningPlasticity.mfront :
   [ parameters      , material ]
   [ Young Modulus   , 70e9     ];
   [ Poisson Ratio   , 0.34     ];
   [ Yield Stress    , 300e6    ];
   [ Hardening Slope , 10e9     ];

   Boundary conditions:
   - uy = 0 on the lower boundary y = 0
   - ux = 0 on the symmetry plane x = 0
   - uz = 0 on the symmetry plane z = 0
   - uy = 6e-3 * t on the upper boundary y = 0.03

   Time: 50 steps from t = 0 to t = 1

   Element:
   - Family H1
   - Order 1

Run the 3D simulation
~~~~~~~~~~~~~~~~~~~~~

Three executables solve this problem with different linear solvers:

- ``ssna303_3d_mumps`` uses MUMPS in parallel and UMFPack sequentially. It
  requires MFEM built with MUMPS. It also offers the FBar formulation.
- ``ssna303_3d_hypre`` uses the FGMRES solver of hypre with the BoomerAMG
  preconditioner.
- ``ssna303_3d_petsc`` uses PETSc with the configuration file ``rc_ex10p``.
  It requires MFEM built with PETSc.

.. code-block:: bash

   mpirun -n 4 ./ssna303_3d_mumps
   mpirun -n 4 ./ssna303_3d_hypre
   mpirun -n 4 ./ssna303_3d_petsc

Available options
~~~~~~~~~~~~~~~~~

+--------------------------------------+-----------------------------------+------------------------+
| Command line                         | Description                       | Default                |
+======================================+===================================+========================+
| ``--order`` or ``-o``                | Finite element order              | 1                      |
+--------------------------------------+-----------------------------------+------------------------+
| ``--nbsteps`` or ``-ns``             | Number of time steps              | 50                     |
+--------------------------------------+-----------------------------------+------------------------+
| ``--end-time`` or ``-et``            | End time. The displacement of the | 1                      |
|                                      | upper boundary is 6e-3 * t.       |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--reference-file`` or ``-rf``      | Reference values of the resultant | no comparison          |
|                                      | force on the upper boundary       |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--use-fbar`` or ``-fb``,           | Use or not the FBar formulation.  | no FBar                |
| ``--no-use-fbar`` or ``-no-fb``      | It requires MGIS built with TFEL. |                        |
|                                      | Only for ``ssna303_3d_mumps``.    |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--standard-reference-file`` or     | Reference values computed without | no comparison          |
| ``-srf``                             | FBar, compared with a larger      |                        |
|                                      | tolerance. Only for               |                        |
|                                      | ``ssna303_3d_mumps``.             |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--parallel`` or ``-p``,            | Run in parallel with MUMPS or     | parallel               |
| ``--no-parallel`` or ``-no-p``       | sequentially with UMFPack. Only   |                        |
|                                      | for ``ssna303_3d_mumps``.         |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--linearsolver`` or ``-ls``        | Linear solver. Only for           | HypreFGMRES            |
|                                      | ``ssna303_3d_hypre``.             |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--preconditioner`` or ``-pc``      | Preconditioner of the linear      | HypreBoomerAMG         |
|                                      | solver. Only for                  |                        |
|                                      | ``ssna303_3d_hypre``.             |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--refinement`` or ``-r``           | Number of uniform refinements of  | 0                      |
|                                      | the mesh. Not for                 |                        |
|                                      | ``ssna303_3d_mumps``.             |                        |
+--------------------------------------+-----------------------------------+------------------------+
| ``--use-petsc`` and                  | Use PETSc with the given          | ``rc_ex10p`` for       |
| ``--petsc-configuration-file``       | configuration file. It requires   | ``ssna303_3d_petsc``,  |
|                                      | MFEM built with PETSc. Not for    | no PETSc otherwise     |
|                                      | ``ssna303_3d_hypre``.             |                        |
+--------------------------------------+-----------------------------------+------------------------+

Tests
~~~~~

In the full test mode, ``ssna303_3d_mumps`` and ``ssna303_3d_hypre`` compute
the whole loading in 10 time steps. In the restricted test mode, they only
compute the first time step, up to 0.1. The resultant force is compared to
``ssna303_3d-force-10steps.ref``.

The test ``ssna303_3d_mumps-fbar`` does the same with FBar. It compares the
resultant force to ``ssna303_3d-force-10steps-fbar.ref``, and to the reference
values without FBar with a larger tolerance.

The test ``ssna303_3d_petsc`` always computes one time step up to 0.02. It
compares the resultant force to ``ssna303_3d-force-1step.ref``.
