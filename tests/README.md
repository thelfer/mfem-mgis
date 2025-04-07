# A brief description of the unit-tests

## Mechanical tests

- `UniaxialTensileTest`: a simple test on a unit cube to test various
  behaviours:
  - an orthotropic elastic behaviour,
  - a simple isotropic behaviour,
  - one of the Mazars damage behaviours,
  - the Saint-Venant-Kirchhoff hyperelastic behaviour
  The results are compared to reference results
- `RveNonLinearElastic`: a test of a Representative Volume Element using an elastic behavior law
- `RveNonLinearElasticV2`: same test using another constructor of the PeriodicNonLinearEvolutionProblem class
- `ImposedPressureTest`: a simple test to check that the imposed
  pressure boundary condition works as expected.
- `MixedOxideFuels`: a simple test of a RVE using an elasto-viscoplastic law, the full simulation is available at `ex7` in the `mfem-mgis-examples` github repository.
- `SatohTest`: this test models a 2D plate of lenght 1 in plane strain
  clamped on the left and right boundaries and submitted to a parabolic
  thermal gradient along the x-axis:
  - the temperature profile is minimal on the left and right boundaries
  - the temperature profile is maximal for x = 0.5
  This example shows how to define an external state variable using an
  analytical profile.
	   
## Heat transfer

- `StationaryNonLinearHeatTransferTest`: a simple test checking that the
  stationary nonlinear heat transfer behaviour integrator works as
  exepected.

## Micromorphic damage

- `MicromorphicDamage2DTest`: this test compares the calculation make
  with the micromorphic with damage behaviour integrator a manufacturaed
  solution in 2D.
- `MicromorphicDamage2DTest2`: this test set describes the failure of a
  bar using an alternate minimization algorithm between a mechanical
  model and a model describing damage evolution in 2D.
- `MicromorphicDamage3DTest`: this test set describes the failure of a
  bar using an alternate minimization algorithm between a mechanical
  model and a model describing damage evolution in 2D. solution.

## Unit tests
- `PeriodicTest`: This test provides some tests of periodic features 
- `ParallelReadMode`: This test checks that the reader can read correctly a splitted mesh
- `PartialQuadratureFunctionTest`:  This test provides some tests
  on partial quadrature functions.



