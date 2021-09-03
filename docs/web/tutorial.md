---
title: A tutorial introduction to `MFEM/MGIS`
author: Thomas Helfer, Guillaume Latu
date: 2021
lang: en-EN
link-citations: true
colorlinks: true
figPrefixTemplate: "$$i$$"
tabPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
---

This tutorial describes how to describe in `MFEM/MGIS` a tensile test on
a notched beam made of an isotropic plastic behaviour with linear
hardening in the logarithmic space. This tutorial highlights the key
features of this project.

The full code is available in the `mfem-mgis-examples` repository in `ex3` directory:

<https://github.com/latug0/mfem-mgis-examples>

# Description of the test case

This tutorial considers a tensile test on a notched beam which is
modelled by plastic behaviour at finite strain:

- The geometry and the mesh is described in Section
  @sec:mfem_mgis:ssna303:mesh.
- The boundary conditions are described in Section
  @sec:mfem_mgis:ssna303:bc.
- The mechanical behaviour is described in Section
  @sec:mfem_mgis:ssna303:behaviour.

This case is a variant of another one available on `code-aster`web site:

<https://www.code-aster.org/V2/doc/v10/fr/man_v/v6/v6.01.303.pdf>

## Geometry and mesh {#sec:mfem_mgis:ssna303:mesh}

![Mesh used to described the notched beam](img/ssna303-mesh.png){#fig:mfem_mgis:ssna303:mesh width=30%}

For symmetry reasons, only half of the notched beam is represented in
Figure @fig:mfem_mgis:ssna303:mesh. The height \(h\) of the beam is
\(30\,mm\). The half-width \(w\) of the beam is \(5.4\,mm\).

The positions of the points \(p1\) \(p2\) and \(c\) are respectivly
\(\left(3\,mm, 0\right)\), \(\left(5.4\,mm, 4.8\,mm\right)\) and
\(\left(9\,mm, 0\right)\).

This notched beam has been meshed using
[`̀Cast3M`](http://www-cast3m.cea.fr/) and exported in the `MED` file
format proposed and used by [Salomé platform](https://www.salome-platform.org/).
This file has been converted in the `msh` file format using
[`gmsh`](https://gmsh.info/) tool in order to import easily in MFEM.

> **Support of the `MED` file format**
>
> Direct support of the `MED` file format is currently under
> development.

## Modelling hypothesis

The beam is treated using the plane strain modelling hypothesis. In
finite strain, this assumes that the axial component of the deformation
gradient is set equal to \(1\).

## Boundary conditions {#sec:mfem_mgis:ssna303:bc}

Dirichlet boundary conditions force the solution to attain certain a priori prescribed values on some boundaries.
The beam is fixed along the bottom line (\(y=0\)\) and a displacement
\(U_{y}\) is imposed an the top of the beam (\(y=h\)).

The symmetry axis on the left is blocked in the `x`-direction.

## Mechanical behaviour {#sec:mfem_mgis:ssna303:behaviour}

### Description

The material of the notched beam is described by a simple isotropic
elasto-plastic behaviour with isotropic hardening in the logarithmic
space [@miehe_anisotropic_2002] and is implemented using the [`MFront`
code generator](http://tfel.sourceforge.net).

This behaviour is characterized by four parameters:

- The `Young Modulus` ($E$) is the slope of the linear part of the
  stress-strain curve for a material under tension or compression
  (isotropic elastic material).
- The `Poisson Ratio` ($\nu$) is the coefficient to characterize the
  contraction of the material perpendicular to the direction of the
  force applied.
- The `Yield Strength` ($\sigma_{0}$) defines the point on the
  stress versus strain curve where the material initially starts to go
  into plastic strain.
- The `Strain Hardening Modulus` (H) defines the slope of the
  stress versus strain curve after the point of yield of a material.

In our example the following values are used:

\[
\left\{
    \begin{array}{lcl}
        \nu & = & 0.34 \\
        \epsilon & = & 70.10^{9} MPa \\
        H & = & 10.10^{9} \\
        s_0 & = & 300.10^{6}
    \end{array}
\right.
\]


### Compilation of the `MFront` behaviour

The previous values are hard-coded in the `MFront` file.
The `MFront` implementation is stored in a source file called
`Plasticity.mfront`. This file must be compiled before the execution of
our `MFEM/MGIS` `C++` example which will be detailed in depth in
Section @sec:mfem_mgis:ssna303. Compilation is performed as follows:

~~~~{.cxx}
$ mfront --obuild --interface=generic Plasticity.mfront
Treating target : all
The following library has been built :
- libBehaviour.so :  Plasticity_AxisymmetricalGeneralisedPlaneStrain 
  Plasticity_Axisymmetrical Plasticity_PlaneStrain
  Plasticity_GeneralisedPlaneStrain Plasticity_Tridimensional
~~~~

# Numerical resolution {#sec:mfem_mgis:ssna303}

## Initialization of the resolution

The `initialize` function must be called at the very beginning of the
`main` function to process the command line arguments:

~~~~{.cxx}
mfem_mgis::initialize(argc, argv);
~~~~

> **The `mfem_mgis` namespace**
>
> All the classes and funtions of the `MFEM/MGIS` project are place in
> the `mfem_mgis` namespace.

This call is mostly useful in parallel and handles:

- The initialization of interprocess communications handled by the
  [`MPI` framework](https://www.mpi-forum.org/docs/).
- The initialization of the [`PETSc` scientific
  toolkit](https://www.mcs.anl.gov/petsc/), if supported and requested.

### Constant variables

The code then defines some constant variables defining the path to the
mesh file, the path the the `MFront` shared library, and the name of the
behaviour:

~~~~{.cxx}
  const char* mesh_file = "ssna303.msh";
  const char* library = "src/libBehaviour.so";
  const char* behaviour = "Plasticity";
~~~~

### Command line options

The numerical resolution can be parametrized using command line
options by relying on the `MFEM` facilities provided by the
`OptionsParser` class.

The proposed implementation allows the following options:

- `--order` which specifies the finite element order (polynomial degree).
- `--parallel` which specifices if the simulation must be run in parallel.

Those options options are associated with local variables which are
default initialized as follows:

~~~~{.cxx}
  auto order = 1;
#if defined(MFEM_USE_MPI)
  bool parallel = true;
#else
  bool parallel = false;
#endif
~~~~

If left unchanged, those default values select:

- a parallel computation if `MFEM` was built with `MPI` support and a
  sequential computation otherwise.
- the use of linear elements.

If `MFEM` was built with support of `PETSc` library, the following
options are added by the `mfem_mgis::declareDefaultOptions` function:

- `--use-petsc` which speficies that linear and non linear solvers of
  the `PETSc` toolkit must be used.
- `--petsc-configuration-file` which specifies a configuration file for
  the `PETSc` toolkit.

In practice, an object of the class `mfem::OptionsParser` is declared.
The expected options are declared and the `Parse` method is called:

~~~~{.cxx}
  mfem::OptionsParser args(argc, argv);
  mfem_mgis::declareDefaultOptions(args);
  args.AddOption(&parallel, "-p", "--parallel",
                 "Perform parallel computations.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return EXIT_FAILURE;
  }
~~~~

## Declaring the non linear problem

The non linear evolution problem is defined as follows:

~~~~{.cxx}
 mfem_mgis::NonLinearEvolutionProblem problem(
      {{"MeshFileName", mesh_file},
       {"FiniteElementFamily", "H1"},
       {"FiniteElementOrder", order},
       {"UnknownsSize", dim},
       {"Hypothesis", "PlaneStrain"},
       {"Parallel", parallel}});
~~~~

The constructor of the `NonLinearEvolutionProblem` class takes an object
of `Parameters` type which is able to store various kind of data in a
hierarchical structure. The valid parameters for the construction of a
non linear evolution problem are described in the `doxygen`
documentation of the `NonLinearEvolutionProblem` class.

The `NonLinearEvolutionProblem` class is the main class manipulated by
the end-users of the `MFEM/MGIS` library. It is meant to handle all the
aspects of the non linear resolution.

Thanks to the `Parameters` type, which is used at different locations in the
interface of the `NonLinearEvolutionProblem` class, the `MFEM/MGIS`
exposes a high level API (Application Programming Interface) which hides
(by default) all the details related to parallelization and memory
management. For example, the parameter `Parallel` allows to switch from
a parallel computation to a parallel one at runtime.

> **Input files and `python` wrappers**
> 
> This high level API can be used to configure a resolution from an
> input file or to wrap the library in `python`. Those features are not
> yet implemented.

Although based on the `MFEM` library, the standard end-user of the
`MFEM/MGIS` library would barely never used directly the `MFEM`
data-structures. However, the `MFEM/MGIS` library does not preclude
to directly use the `MFEM` data-structures, built-in
non linear forms, etc. This lower level API is not described in this
tutorial.

## Names boundaries and materials

`MFEM` distinguishes elements of the mesh (materials and boundaries) by
integers. This may seem unpractical to most users. The `MFEM/MGIS` allows
to associate names to materials and boundaries as follows:

~~~~{.cxx}
  problem.setMaterialsNames({{1, "NotchedBeam"}});
  problem.setBoundariesNames(
      {{3, "LowerBoundary"}, {4, "SymmetryAxis"}, {2, "UpperBoundary"}});
~~~~

> **Automatic definition of the names of materials and boundaries**
>
> Many mesh file formats naturally associate names to mesh elements.
> This is the case for `MED` file format and the `msh` file format
> generated by `gmsh`.
>
> Future versions of the library may thus automatically define the names
> of materials and boundaries.

## Declaring the mechanical behaviour

The following line associates a mechanical behaviour to the first
material:

~~~~{.cxx}
  problem.addBehaviourIntegrator("Mechanics", "NotchedBeam",
                                 "src/libBehaviour.so",
`                                 "Plasticity");
~~~~

The four arguments of the `addBehaviourIntegrator` are:

- The type of physical problem described. Currently two types of
  physical problems are supported out of the box by the library:
  `Mechanics` and `HeatTransfer`. Support for other physical problems
  can be plugged in at runtime if needed.
- The material identifier, as defined in the mesh file. This identifier
  may be either an integer or a string. In the later case, the string is
  interpreted as a regular expression, a feature introduced by the
  `Licos` fuel performance code and which proved very pratical in many
  cases [@helfer_licos_2015].
- The shared library containing the behaviour to be used.
- The name of the behaviour to be used.

> **Information associated with the behaviour and automatic memory management**
>
> Thanks to the [`MGIS`
> project](https://thelfer.github.io/mgis/web/index.html)
> [@helfer_mfrontgenericinterfacesupport_2020], all the information
> related to the mechanical behaviour is retrieved, including:
>
> - The type of behaviour (finite strain mechanical behaviour is this case).
> - The names of material properties, parameters, state variables and
>   external state variables.
> - etc.
>
> The memory required to store the state of the materials is
> automatically allocated.

## Initialisation of the temperature

The following lines define an uniform temperature on the material at the
beginning of the time step and at the end of time step:

~~~~{.cxx}
  auto& m1 = problem.getMaterial("NotchedBeam");
  mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
  mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
~~~~

Defining the temperature is required by all `MFront` behaviours.

The object returned a by the `getMaterial` method returns a thin wrapper
around the `MaterialDataManager` provided by the [`MGIS`
project](https://thelfer.github.io/mgis/web/index.html)
[@helfer_mfrontgenericinterfacesupport_2020].

In the previous lines, `m1.s0` and `m1.s1` denotes respectively the
state of the material at the beginning of the time step and at the end
of the time step.

## Boundary Condition

The `NonLinearEvolutionProblem` class allows to define uniform Dirichlet
boundary conditions (imposed displacement) using the
`addUniformDirichletBoundaryCondition` method as follows:

~~~~{.cxx}
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "LowerBoundary"}, {"Component", 1}});
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "SymmetryAxis"}, {"Component", 0}});
  problem.addUniformDirichletBoundaryCondition(
      {{"Boundary", "UpperBoundary"},
       {"Component", 1},
       {"LoadingEvolution", [](const auto t) {
          const auto u = 6e-3 * t;
          return u;
        }}});
~~~~

Again, the code is almost self-explanatory. If the value of the imposed
displacement is not specified (using the `LoadingEvolution` parameter),
the selected component is set to zero. The `LoadingEvolution` parameter
allows to specify the evolution of the imposed displacement using a
function of time (defined her using a `C++` lambda expression).

## Non linear solver parameters.

If `PETSc` is not used, the following line set the parameters of the
Newton-Raphson solver used to find the equilibrium of the whole
structure:

~~~~{.cxx}
  if (!mfem_mgis::usePETSc()) {
    problem.setSolverParameters({{"VerbosityLevel", 0},
                                 {"RelativeTolerance", 1e-6},
                                 {"AbsoluteTolerance", 0.},
                                 {"MaximumNumberOfIterations", 10}});
  }
~~~~

Valid parameters for the `setSolverParameters` are described in the
`doxygen` documentation of the library.

If `PETSc` is used (see the `--use-petcs` command line option), the
parameters associated with the choice of the non linear solver must be
provided by an external configuration file (see the
`--petsc-configuration-file` command line option).

## Selection of the linear solver

If `PETSc` is not used, the linear solver can be selected using the
`setLinearSolver` method. Here we select `MUMPS`, in parallel and
`UMFPack` in sequential:

~~~~{.cxx}
  if (!mfem_mgis::usePETSc()) {
    if (parallel) {
      problem.setLinearSolver("MUMPSSolver", {});
    } else {
      problem.setLinearSolver("UMFPackSolver", {});
    }
  }
~~~~

The second argument is an object of the `Parameters` type which can be
used to fine tune the linear solver and, in the case of iterative
solvers, optionnaly define a preconditioner. For direct solvers, no
parameters are required.

## Post-processings

The `addPostProcessing` method let the user define some built-in
postprocessings.

In this example, we export the displacements for visualization in
[`paraview`](https://www.paraview.org/) and compute the resultant force
on the boundary where the displacement as follows:

~~~~{.cxx}
  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "ssna303-displacements"}});
  problem.addPostProcessing("ParaviewExportIntegrationPointResultsAtNodes",
                            {{{"Results", "FirstPiolaKirchhoffStress"},
                              {"OutputFileName", "ssna303-stress"}}});
  problem.addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{{"Results", "EquivalentPlasticStrain"},
        {"OutputFileName", "ssna303-equivalent-plastic-strain"}}});
  problem.addPostProcessing("ComputeResultantForceOnBoundary",
                            {{"Boundary", 2}, {"OutputFileName", "force.txt"}});
~~~~

These post-processings are called using the `executePostProcessings`
method during the runtime using the state at the end of time step.
The user may also plugged in its own post-processing.

## Resolution

The `NonLinearEvolutionProblem` class is meant to solve the problem on
one time step only. This allows to easily built weakly coupled non
linear resolutions (for example, thermo-mechanical resolutions where the
heat tranfer and mechanical problems are solved using a staggered
scheme) or set-up couplings with external solvers.

In this tutorial, a local time-substepping scheme is set up to handle
resolution failures.

The loading starts at time \(0\) and ends at time \(1\). This range is
divided in \(50\) time steps.

~~~~{.cxx}
  const auto nsteps = mfem_mgis::size_type{50};
  const auto dt = mfem_mgis::real{1} / nsteps;
  auto t = mfem_mgis::real{0};
  auto iteration = mfem_mgis::size_type{};
  for (mfem_mgis::size_type i = 0; i != nsteps; ++i) {
    std::cout << "iteration " << iteration << " from " << t << " to " << t + dt
              << '\n';
~~~~

The local time substepping scheme is simply set up as follows:

~~~~{.cxx}
    auto ct = t;
    auto dt2 = dt;
    auto nsteps = mfem_mgis::size_type{1};
    auto nsubsteps  = mfem_mgis::size_type{0};
    while (nsteps != 0) {
      auto converged = problem.solve(ct, dt2);
      if (converged) {
        --nsteps;
        ct += dt2;
        problem.update();
      } else {
        nsteps *= 2;
        dt2 /= 2;
        ++nsubsteps;
        problem.revert();
        if (nsubsteps == 10) {
          mfem_mgis::raise("maximum number of substeps");
        }
      }
    }
~~~~

Every time a resolution is sucessful, the material state is updated
using the `update` method, the current time is incremented and the
number of the remaining substeps is decreased. The loop stops when the
remaining number of sub-steps goes to zero.

If the resolution failed, the local time step is divided by \(2\), the
number of remaining substeps is multiplied by \(2\) and the state of the
material is reverted to the beginning of the time step using the
`revert` method. The resolution stops if more than \(10\) nested reverts
are generated.

Once a time step has been successful, the post-processings are executed
and the time is incremented.

~~~~{.cxx}
      problem.executePostProcessings(t, dt);
      t += dt;
      ++iteration;
    }
  }
~~~~  

