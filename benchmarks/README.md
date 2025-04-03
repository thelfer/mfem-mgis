# Benchmarks section


Please, refer to the docuementation website for more information: https://thelfer.github.io/mfem-mgis/developer_guide/benchmark.html

It includes:

- Description of the profiling tools
- Benchmarks:
  - RVE of a MOX simulation from example 7 of `mfem-mgis-examples`: https://github.com/latug0/mfem-mgis-examples/tree/master/ex7
  - Overhead between MFEM and MFEM-MGIS 

## Overhead between MFEM and MFEM-MGIS

This benchmark aims to evaluate the time overhead between `MFEM` and `MFEM/MGIS` for a simple test case involving an elastic behavior law.

Details:

```
    Solver: HyprePCG

    Preconditionner: HypreDiagScale (Jacobi)

    Tolerance: 1e-14

    Boundary conditions: we impose an uniform dirichlet condition U = (0,0,0) at left and U = (0,0,1) at right.
```

Note that due to the small size of the input mesh, domain decomposition is limited to 16 subdomains (MFEM-MGIS only performs parallel refinement). To improve this benchmark, consider using sequential refinement or a finer mesh.

Run this example on Topaze supercomputer:

```
ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5
ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMMGISLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5
```
