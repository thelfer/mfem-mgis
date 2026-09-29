Mechanical tests
================

- `UniaxialTensileTest`: a simple test on a unit cube to test various
  behaviours:
  - an orthotropic elastic behaviour,
  - a simple isotropic behaviour,
  - one of the Mazars damage behaviours,
  - the Saint-Venant-Kirchhoff hyperelastic behaviour
  The results are compared to reference results
- `ImposedPressureTest`: a simple test to check that the imposed
  pressure boundary condition works as expected.
	   
Heat transfer
=============

- `StationaryNonLinearHeatTransferTest`: a simple test checking that the
  stationary nonlinear heat transfer behaviour integrator works as
  exepected.

Micromorphic damage
===================

- `MicromorphicDamage2DTest`: this test compares the calculation make
  with the micromorphic with damage behaviour integrator a manufacturaed
  solution in 2D.
- `MicromorphicDamage2DTest2`: this test set describes the failure of a
  bar using an alternate minimization algorithm between a mechanical
  model and a model describing damage evolution in 2D.
- `MicromorphicDamage3DTest`: this test set describes the failure of a
  bar using an alternate minimization algorithm between a mechanical
  model and a model describing damage evolution in 2D. solution.

Point-wise models
=================

- The `PointWiseModelTest`: this test shows how to evaluate a simple
  model of shrinkage of Uranium dioxide under irradiation on a partial
  quadrature space

Unit tests
==========

- `PostProSubMesh`: This test provides some tests of ParaviewExportResults 
- `PeriodicTest`: This test provides some tests of periodic features 
- `ParallelReadMode`: This test checks that the reader can read
  correctly a splitted mesh.
- `PartialQuadratureSpaceTest`: This test provides some tests on
  partial quadrature spaces.
- `PartialQuadratureFunctionTest`: This test provides some tests on
  partial quadrature functions.
- `GridFunctionTest`: test the `update` function which projects a
  `GridFunction` on a partial quadrature function.
- `L2ProjectionTest`: test the `computeL2Projection` function which
  projects a partial quadrature function on a `GridFunction`.
- `ImplicitGradientRegularizationTest`: test the
  `computeImplicitGradientRegularization` which computes the implicit
  gradient regularization of a partial quadrature function.

<!--
## Benchmark

- `elasticity`: pure elastic test base 

-->