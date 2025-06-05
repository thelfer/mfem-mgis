===========================
Mechanical post-processings
===========================

The :cxx:`MechanicalPostProcessings.hxx` header declares various utility
functions:

- :cxx:`computeVonMisesEquivalentStress`: computes the von Mises
  equivalent stress. For finite strain behaviours, the von Mises
  equivalent stress of the Cauchy stress is returned.
- :cxx:`computeEigenStresses` computes the eigen values of the stress.
  For finite strain behaviours, the eigen values of the Cauchy stress is
  returned.
- :cxx:`computeFirstEigenStress` computes the first (maximum) eigen
  value of the stress. For finite strain behaviours, the first eigen
  value of the Cauchy stress is returned.
- :cxx:`computeStressInGlobalFrame`: compute the stress in the global
  frame. This function is only for orthotropic behaviours. For finite
  strain behaviours, the first Piola-Kirchhoff stress is returned.
- :cxx:`computeCauchyStressInGlobalFrame`: compute the Cauchy stress in
  the global frame. This function is only valid for orthotropic finite
  strain behaviours.
