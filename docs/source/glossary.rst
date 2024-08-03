======================
Glossary Of Parameters
======================

``MFEM-MGIS`` offers a simplified API based on encapsulating ``Parameter`` objects in a ``Parameters`` class, which is then given to various functions to implement the problem to be solved.
It is important to note that if a Parameter is not listed in the glossary, the code will stop and an error message will specify the available options.
In this section we describe the options available for each part of the code, both for setting the parameters of the linear solver and for defining the problem to be solved.

.. note:: 

  The list of default parameters is adapted according to compilation flags. For example, if ``MFEM_USE_MUMPS`` is not active, ``MUMPS`` will not appear in the list. So if one of the options listed below doesn't appear,  ``MFEM-MGIS`` may not have been compiled with all the options.

Setup Your Linear Solver
========================

In this section, we'll distinguish between two types of solver: direct solvers, which are not always listed in the parallel version of the problems, and Krylov solvers, which are available in both serial and parallel versions (except for ``hypre`` solvers). 

To choose your solver, pass a Parameters object to the `setLinearSolver` function.

Direct Solver
-------------

MUMPS Solver
^^^^^^^^^^^^

- Solver name: ``MUMPSSolver``
- Dependencies: 
  - ``MFEM_USE_MUMPS``
  - parallel
- Parameters:
  - ``Symetric`` [bool]: specify if the matrix is symetric.
  - ``PositiveDefinite`` [bool]: specify if the matrix is positive definite

**Example:**

.. code-block:: cpp
  
  problem.setLinearSolver("MUMPSSolver", {{"Symmetric", true}})

UMFPack Solver
^^^^^^^^^^^^^^

- Solver name: ``UMFPackSolver``
- Dependencies: 
  - ``MFEM_USE_SUITESPARSE``
  - not parallel
- Parameters:None

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("UMFPackSolver", {})

.. note::

  Direct solvers do not need any tolerance parameter or preconditioner parameter.

.. note::

  We advise you not to use the direct solvers using the parallel option, this advice is is mainly due to the fact that these solvers don't scale in memory on large problems (e.g. meshes with more than a hundred thousand points).  Nevertheless, these solvers are maintained within the framework of small computations achievable on a standard laptop.

Krylov Solver
-------------

In this section, we'll start with the solvers present in ``mfem`` and then describe the ``hypre`` solvers available only in the parallel version.

Conjugate Gradient
^^^^^^^^^^^^^^^^^^

- Solver name: ``CGSolver``
- Dependencies: None
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner.
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``RelativeTolerance`` [double]: Define the relative tolerance
  - ``AbsoluteTolerance`` [double]: Define the absolute tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("CGSolver", {{"VerbosityLevel", 1},
            {"AbsoluteTolerance", 1e-12},
            {"RelativeTolerance", 1e-12},
            {"MaximumNumberOfIterations", 5000}
            })


Generalized Minimal Residual (GMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``GMRESSolver``
- Dependencies: None
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner.
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``RelativeTolerance`` [double]: Define the relative tolerance
  - ``AbsoluteTolerance`` [double]: Define the absolute tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("GMRESSolver",
          {{"VerbosityLevel", 1},
          {"AbsoluteTolerance", 1e-12},
          {"RelativeTolerance", 1e-12},
          {"MaximumNumberOfIterations", 100000}});


Biconjugate Gradient Stabilized (BiCGSTAB)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``BiCGSTABSolver``
- Dependencies: None
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner.
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``RelativeTolerance`` [double]: Define the relative tolerance
  - ``AbsoluteTolerance`` [double]: Define the absolute tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("BiCGSTABSolver",
          {{"VerbosityLevel", 1},
          {"AbsoluteTolerance", 1e-12},
          {"RelativeTolerance", 1e-12},
          {"MaximumNumberOfIterations", 1000}});


Minimal Residual (MINRES)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``MINRESSolver``
- Dependencies: None
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner.
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``RelativeTolerance`` [double]: Define the relative tolerance
  - ``AbsoluteTolerance`` [double]: Define the absolute tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("MINRESSolver",
          {{"VerbosityLevel", 1},
          {"AbsoluteTolerance", 1e-12},
          {"RelativeTolerance", 1e-12},
          {"MaximumNumberOfIterations", 1000}});

Stationary Linear Iteration (SLI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``SLISolver``
- Dependencies: None
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner.
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``RelativeTolerance`` [double]: Define the relative tolerance
  - ``AbsoluteTolerance`` [double]: Define the absolute tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example:**

.. code-block:: cpp

  problem.setLinearSolver("SLISolver",
          {{"VerbosityLevel", 1},
          {"AbsoluteTolerance", 1e-12},
          {"RelativeTolerance", 1e-12},
          {"MaximumNumberOfIterations", 1000}});

Preconditioned Conjugate Gradient (HyprePCG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HyprePCG``
- Dependencies:
  - ``hypre``
  - parallel
- Parameters:
  - ``Precondition`` [Parameters]: Define your preconditioner. 
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``Tolerance`` [double]: Define the tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example with a Jacobi preconditioner:**

.. code-block:: cpp

    problem.setLinearSolver("HyprePCG", {{"VerbosityLevel", 1},
          {"Tolerance", 1e-12},
          {"Preconditioner", diagscale},
          {"MaximumNumberOfIterations", 5000}});

.. note::
  
  Only ``hypre`` preconditioners are allowed with the HyprePCG linear solver.

Generalized Minimal Residual (HypreGMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HypreGMRES``
- Dependencies:
  - ``hypre``
  - parallel
- Parameters:
  - ``KDim`` [int]: Define the Krylov dimension.
  - ``Precondition`` [Parameters]: Define your preconditioner. 
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``Tolerance`` [double]: Define the tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example with a  Parasail preconditioner:**

.. code-block:: cpp

   problem.setLinearSolver("HypreGMRES", {{"VerbosityLevel", 1},
          {"Tolerance", 1e-12},
          {"Preconditioner", parasail},
          {"MaximumNumberOfIterations", 5000}});

.. note::
  
  Only ``hypre`` preconditioners are allowed with the HypreGMRES linear solver.


Flexible GMRES (HypreFGMRES)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Solver name: ``HypreFGMRES``
- Dependencies:
  - ``hypre``
  - parallel
- Parameters:
  - ``KDim`` [int]: Define the Krylov dimension.
  - ``Precondition`` [Parameters]: Define your preconditioner. 
  - ``VerbosityLevel`` [0,1]: Define the verbosity of the solver
  - ``Tolerance`` [double]: Define the tolerance
  - ``MaximumNumberOfIterations`` [int]: Maximum number of Krylov iterations

**Example with a ILU preconditioner:**

.. code-block:: cpp

  problem.setLinearSolver("HypreFGMRES", {{"VerbosityLevel", 1},
          {"Tolerance", 1e-12},
          {"Preconditioner", ilu},
          {"MaximumNumberOfIterations", 5000}});


.. note::
  
  Only ``hypre`` preconditioners are allowed with the HypreFGMRES linear solver.

Setup Your Preconditioner
=========================

Currently, only preconditioners from ``hypre`` are integrated into ``MFEM-MGIS``.

Algebraic MultiGrid (AMG) 
-------------------------

- Preconditioner name: ``HypreBoomerAMG``
- Options:
  - ``Strategy`` [Elasticity, System, None]: Propose strategies that can lead to improve convergence and scalability
  - ``VerbosityLevel`` [0,1]:  Define the verbosity of the preconditioner
- website: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}, {"Strategy", System}};
    auto amg = mfem_mgis::Parameters{{"Name","HypreBoomerAMG"}, {"Options",options}};

**Reference:**

.. code-block:: text

  J. W. Ruge and K. Stüben. Algebraic multigrid (AMG). In S. F. McCormick, editor, Multigrid Methods, volume 3 of Frontiers in Applied Mathematics, pages 73–130. SIAM, Philadelphia, PA, 1987.

Euclid (HypreEuclid)
--------------------
- Preconditioner name: ``HypreEuclid``
- Options:
  - ``VerbosityLevel`` [0,1]:  Define the verbosity of the preconditioner
- website: https://hypre.readthedocs.io/en/latest/solvers-euclid.html 

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto euclid = mfem_mgis::Parameters{{"Name","HypreEuclid"}, {"Options",options}};

**Reference:**

.. code-block:: text

  D. Hysom and A. Pothen. Efficient parallel computation of ILU(k) preconditioners. In Proceedings of Supercomputing ‘99. ACM, November 1999. Published on CDROM, ISBN #1-58113-091-0, ACM Order #415990, IEEE Computer Society Press Order # RS00197.

Incomplete LU factorization (HypreILU)
--------------------------------------

- Preconditioner name: ``HypreILU``
- Dependencies:
  - ``MFEM_HYPRE_VERSION`` >= 21900
- Options:
  - ``HypreILULevelOfFill`` [int]:
  - ``VerbosityLevel`` [0,1]:  Define the verbosity of the preconditioner
- website: https://hypre.readthedocs.io/en/latest/solvers-ilu.html

**Example:**

.. code-block:: cpp

    auto ilu = mfem_mgis::Parameters{{"Name","HypreILU"}, {"Options", mfem_mgis::Parameters{{"HypreILULevelOfFill", 1}}}};


Sparse Approximate Inverse (HypreParaSails)
-------------------------------------------

- Description: ParaSails is a parallel implementation of a sparse approximate inverse preconditioner, using a priori sparsity patterns and least-squares (Frobenius norm) minimization.
- Preconditioner name: ``HypreParaSails``
- Options:
  - ``VerbosityLevel`` [0,1]:  Define the verbosity of the preconditioner
- website: https://hypre.readthedocs.io/en/latest/solvers-parasails.html

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto parasail = mfem_mgis::Parameters{{"Name","HypreParaSails"}, {"Options",options}};

**Reference:**

.. code-block:: text

  E. Chow. A priori sparsity patterns for parallel sparse approximate inverse preconditioners. SIAM J. Sci. Comput., 21:1804–1822, 2000.

.. warning::

  ParaSails is not actively supported by the hypre development team. 

Jacobi (HypreDiagScale)
-----------------------

- Preconditioner name: ``HypreDiagScale``
- Options:
  - ``VerbosityLevel`` [0,1]:  Define the verbosity of the preconditioner

**Example:**

.. code-block:: cpp

    auto options = mfem_mgis::Parameters{{"VerbosityLevel", 0}};
    auto diagscale = mfem_mgis::Parameters{{"Name","HypreDiagScale"}, {"Options",options}};

Setup Your NonLinearProblem
===========================

.. warning::

  This section is under construction.

PeriodicNonLinearProblem
------------------------


Material parameters
-------------------

Material Identifier
^^^^^^^^^^^^^^^^^^^

Functions: getMaterialIdentifier or 

- Key: ``Material`` or ``Materials``
- 

**Example:**

Boundary Condition parameters
-----------------------------



Boundary Mutators:
^^^^^^^^^^^^^^^^^^

- Functions: `getBoundariesIdentifiers`
- Key: ``Boundary`` or ``Boundaries``

**Example:**

Dirichlet Condition
-------------------

This is how to apply a dirichlet boundary condition to several boundaries.

**Example:**

.. code-block:: cpp

  // boundary conditions
  for (const auto boundary : {"left", "right"}) {
    for (const auto dof : {0, 1}) {
      problem.addBoundaryCondition(
          std::make_unique<mfem_mgis::UniformDirichletBoundaryCondition>(
              problem.getFiniteElementDiscretizationPointer(), boundary, dof));
    }
  }

