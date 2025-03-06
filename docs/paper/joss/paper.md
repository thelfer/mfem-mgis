---
title: 'MFEM/MGIS, a HPC mini-application targeting nonlinear thermo-mechanical simulations of nuclear fuels at mesoscale'
tags:
  - thermo-mechanical
  - HPC
  - mfem
  - nuclear fuel
authors:
  - name: Thomas Helfer
    orcid: 0000-0003-2460-5816
    affiliation: 1
  - name: Guillaume Latu
    orcid: 0009-0001-7274-1305
    affiliation: 1
  - name: Raphaël Prat
    orcid: 0009-0002-3808-5401
    affiliation: 1
  - name: Maxence Wangermez
    orcid: 0000-0002-3431-5081
    affiliation: 1
  - name: Francesca Cuteri
    orcid: 0000-0002-9339-8892
    affiliation: 1
affiliations:
 - name: CEA, DES, IRESNE, DEC, Cadarache F 13108 St Paul Lez Durance
   index: 1
date: 6 March 2025
bibliography: paper.bib
---

# Summary

The `MFEM/MGIS` application aims to efficiently use supercomputers in
order to describe coupled multiphysics phenomena with a particular focus
on thermo-mechanics. The authors primarily aim at describing the
nuclear fuels at mesoscale (see example below), but `MFEM/MGIS` is
versatile and can address more general cases.

This open-source library is based on several components as
prerequisites: the `MFEM` (Modular Finite Element Methods) [@mfem;
@mfem-web] library, the `MGIS` (MFront Generic Interface Support)
[@Helfer2020] library and the `MFront` code generator
[@helfer2015introducing].

Thanks to the features embedded within `MGIS` and `MFront` and thanks to
specific developments, `MFEM/MGIS` adds several mechanical features
compared to a pure MFEM approach. The library tackles some peculiarities
of nonlinear mechanics. In particular, the support of complex
constitutive laws and the management of advanced boundary conditions.

# Context and motivation

## About MFEM

`MFEM`, is a finite element library designed for current
supercomputers and the upcoming exascale supercomputers. It
provides many useful features for carrying out state-of-the art
simulations: support for curvilinear meshes, high order approximation
spaces and different families of finite elements, interfaces to
several types of parallel solvers (including matrix-free ones),
preconditioners, and native support for adaptive non-conforming mesh
refinement (AMR).

Originating from the applied mathematics and parallel computing
communities, `MFEM` offers both performance and a large panel of
advanced mathematical features. In particular, one can easily switch
from one linear solver to another (direct or iterative), essential for the targeted application [@bernaud2024pleiades]:
microstructure and mesoscale modeling for nuclear fuel.

## Statement of need

The solid mechanic examples in `MFEM` are mostly limited to simple
constitutive equations such as elasticity and hyperelasticity without
internal state variables. This is insufficient to address many
engineering studies and in particular complex nuclear fuel
simulations.

The aim of `MFEM/MGIS` project is to combine `MFEM` with the
`MFrontGenericInterfaceSupport` (`MGIS`) project, an open-source `C++`
library handles all the kinds of behavior supported by the open-source
`MFront` code generator.

In the field of nonlinear mechanics, this encompasses arbitrary complex
behaviors that can describe damage, plasticity, and viscoplasticity in both
small or finite strain analyses. Generalized behaviors such as
variational approaches to fracture are supported by `MFEM/MGIS`.

The `MGIS` data structures are used to add support for partial
quadrature functions to `MFEM`, a feature needed to store internal state
variables on each material.

## State of the field

Many open-source thermomechanical solvers allow to handle complex
mechanical behaviours. `code_aster`, `MoFEM`, `CalculiX` are examples
of state of the art solvers which have an interface with `MFront`.

However, those solvers lack many features provided by `MFEM` that the
authors found interesting to explore in the field of solid mechanics
(see the above section for a detailed list). The authors also found
interesting to take a platform designed from the start for high
performance computing and adapt it to engineering needs and evaluate
the resulting performances.

<!--
Strongly coupled thermo-mechanical behaviors and Cosserat media can be
introduced by plugins and will be integrated in future versions of
`MFEM/MGIS`.
-->

# Overview of `MFEM/MGIS` features

## The `NonLinearEvolutionProblem` class

The main class of `MFEM/MGIS` is called `NonLinearEvolutionProblem` and
describes the evolution of the materials of the physical system
of interest over a single time step for a given phenomenon.

Currently `MFEM/MGIS` provides built-in support for mechanics, heat
transfer, and micromorphic damage are supported.

The following snippet declares a new nonlinear evolution problem:

~~~~.c++
mfem_mgis::NonLinearEvolutionProblem problem(
    {{"MeshFileName", "fuel.msh"},
     {"FiniteElementFamily", "H1"},
     {"FiniteElementOrder", 6},
     {"UnknownsSize", 1},
     {"Hypothesis", "Tridimensional"},
     {"Parallel", true}});
~~~~

As the unknown is scalar (according to the `UnknownsSize` parameter),
this problem can be used to describe heat transfer or micromorphic
damage, depending on the behaviour integrators declared, as explained in
the next section. The `NonLinearEvolutionProblem` class supports both
sequential and parallel computations and let the user exploits a large
subset of `MFEM` abilities, including the use of finite elements of
arbitrary orders.

A staggered approach for multiphysics simulations can be set up by using
several instances of `NonLinearEvolutionProblem`.

### The `PeriodicNonLinearEvolutionProblem` class

`MFEM/MGIS` provides a specialized version for the
`NonLinearEvolutionProblem` for periodic computations named
`PeriodicNonLinearEvolutionProblem`.

This class allows to manage the evolutions of the macroscopic
gradients (strain in small strain analysis, deformation gradient in
finite strain analysis, temperature gradient in heat transfer
analysis) and pass them to the behaviour integrators.

### Note about linear analyses
 
As implied by its name, the `NonLinearEvolutionProblem` is focused on
nonlinear resolutions. Linear analyses can still be performed by using
linear behaviours (generated by `MFront`), but with a computational
overhead compared to linear analysis made with optimised kernels, such
as the elastic kernels provided natively by `MFEM`.

In our experience, this overhead is limited and mostly comes from
the extra flexibility allowed by `MFEM/MGIS`. For instance,
`MFEM` elastic kernels assume that:

- all materials are elastic,
- material properties (Young's modulus, Poisson's ratio) are uniform
  on each material.

`MFEM/MGIS` kernels don't have those restrictions as described below.

## Behaviour integrators

Contrary to non linear forms provided natively by `MFEM`, `MFEM/MGIS`
allows to assign distinct behaviors to each material. To achieve this, a
special nonlinear formulation has been implemented which delegates the
computations of residual and jacobian terms on each material to
so-called behaviour integrators.

Behaviour integrators are associated with a physical phenomenon and a
modelling hypothesis (plane strain, plane stress).

This following snippet assigns a behaviour integrator to the material
named `beam` to a mechanical non linear evolution problem:

~~~~.c++
mechanics.addBehaviourIntegrator("Mechanics", "beam",
                                  "src/libBehaviour.so",
                                  "MicromorphicDamageI_SpectralSplit");
~~~~

The behaviour `MicromorphicDamageI_SpectralSplit` is loaded from a
behaviour named `libBehaviour.so` which shall have been generated using
`MFront` before running the simulation. The behaviour integrator
`Mechanics` supports arbitrary small strain and finite strain behaviours.

Internally, the `addBehaviourIntegrator` method calls an abstract
factory which instanciates a `BehaviourIntegrator` dedicated to the kind
of behaviour selected by the user (small or finite strain) and the
modelling hypothesis declared by the problem (plane strain, plane
stress, tridimensional, etc.).

### About the definition of material properties

Behaviours may require the user to provide properties, such as the
Young's modulus, Poisson's ratio, etc.. In `MFEM/MGIS`, those
properties can be uniform on the material or given by a partial
quadrature function. The later case allows the properties to be
defined independently on each integration point, which is required if
those properties depend on local material properties or the local
state of the material (for instance, the local temperature).

## User interface

The `MFEM/MGIS` library is written in `C++17` language.

As the application targets mechanical engineers, it provides a high level
of abstraction, focused on the physical aspects of the simulation and
hiding most numerical details by default.

The API is declarative and mostly based on data structures similar to a
`python` dictionary, limiting direct usage of `C++`. In particular, this
data structure is used to instantiate non linear evolution problems,
behaviour integrators, post-processings, and boundary conditions.

<!--
This data structure can be read from a `JSON`-file allowing to create
domain-specific applications.
-->

## Post-processings

Various post-processings are available. Here are some examples of
post-processings that were added to `MFEM/MGIS`:

- `ComputeResultantForceOnBoundary`: Compute the resultant of the inner
  forces on a boundary.
- `ComputeMeanThermodynamicForcesValues`: Compute the macroscopic stress
  and strain for each material. This is mostly useful for 
- `ParaviewExportIntegrationPointResultsAtNodes`: Paraview post
  processing files of partial quadrature functions, like the ones
  associated with the internal state variables.

## Tutorials

These features are described in the following tutorial:
<https://thelfer.github.io/mfem-mgis/tutorial.html>.

Several examples can be found on the open-source GitHub repository:
<https://github.com/latug0/mfem-mgis-examples>, including simulation of
the microstructural scale using a Representative Volume Element (RVE) of
nuclear fuel, figure \autoref{fig:RVE}. 


![Simulation of a Representative Volume Element (RVE) of mixed oxide material, which contains 634 inclusions. The variable "u" represents the displacement. \label{fig:RVE}](./MOX17percent634Spheres.png){width=70%}


<!--
## Plugins

Post-processings and boundary conditions are instantiated through
abstract factories. Those abstract factories allow to extend `MFEM/MGIS`
through plugins.
-->

## Software stack and installation process

![MFEM/MGIS Stack. Each library is open source.\label{fig:SoftStack}](./MMM-stack-v0.png "MFEM/MGIS software stack. MFEM/MGIS is the convergence between two open sources ecosytem."){width=75%}

The MFEM/MGIS software stack is illustrated in \autoref{fig:SoftStack}.
Hence, the minimal package requirements to build MFEM/MGIS on
a HPC platform is typically:`C++17`, `MFEM`, `MGIS`, `TFEL`(MFront),
CMake and `MPI`. It is important to note that `MFEM` includes a significant stack of several tens packages including linear solver libraries.
To handle the numerous accessible combinations, `Spack`
[@gamblin2015spack] is really a cornerstone. This package manager
simplifies building, installing, customizing, and sharing HPC software
stacks. In the end, `MFEM-MGIS` offers a `Spack`
package reusing most of the `MFEM` installation variants already provided in
the `MFEM` `Spack` package while maintaining compatibility across package versions.


<!--

# Integrate MFEM/MGIS in an Open Source Ecosystem

The `MFEM/MGIS` library takes full advantage of an open-source software
(OSS) stack. It benefits from the increasing maturity of several
scientific tools that are combined.

Thus, within `MFEM`, one has many available choices to set the linear
solver, such as: `Hypre`[@hypre2002], `PETSc`[@petsc-web-page],
`MUMPS`[@mumps], `SuperLU`, `UMFPACK`[@davis2004algorithm] or other
ones. Likewise, several preconditioners, partitioning libraries, or
input mesh formats can be activated and used.

Combinations are highly configurable and almost all external libraries
that relates to linear solvers are switchable.
-->

# Numerical Results {#sec:numerical_results}

 Multi-material elastic modelling on computational clusters has been
carried out with `MFEM/MGIS`. The observed scalability performance is
good on a few thousands of CPU cores. Benchmarks are available in the
documentation: <https://thelfer.github.io/mfem-mgis/benchmark.html>.

Several examples can be found on the open-source GitHub repository:
<https://github.com/latug0/mfem-mgis-examples>.

<!--

# Conclusion

This paper presents the `MFEM/MGIS` HPC application designed to address
large scale thermo-mechanical simulation and recent supercomputers.
Based on an open source software stack, it allows for the fine
representation of microstructure in full 3D in the field of fuel
modeling.

On the one hand, `MGIS` and `MFront` bring support for complex
nonlinear behaviours such as damage, plasticity, viscoplasticity
capabilities. On the other hand, `MFEM` provides advanced finite
element schemes and parallel performance (tested on several thousands
of cores until now). 
-->


<!--
Regarding performance portability on GPUs, MFEM already offers numerous
algorithms such as partial assembly on GPUs, but MFEM/MGIS does not
exploit these features yet. Work is underway to port behavior laws to
the GPU (`MFront`) and associated data structures (`MGIS`). Note that
other HPC libraries is considered for the future, such as `CUDA`, `RAJA`
or `OCCA` for performance portability of MFEM/MGIS on GPU.

-->

# Acknowledgement

Funded by the European Union, this work has received funding from the
Open HPC thermomechanical tools for the development of EATF fuels
undertaking (OperaHPC) grant agreement No 101061453.

Benchmarks and scalability tests were performed using HPC resources from
CCRT funded by the CEA/DES simulation program.

# References


