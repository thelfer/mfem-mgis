#! /usr/bin/env bash

# compilation of the code generator
g++ -Wall -W -pedantic -Werror -std=c++17 -g behaviour-integrator.cxx -o behaviour-integrator -lginac -lcln
#
./behaviour-integrator -s isotropic --hypothesis=Tridimensional --generator=StandardFiniteStrainMechanics --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=Tridimensional --generator=StandardSmallStrainMechanics --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=PlaneStrain --generator=StandardFiniteStrainMechanics --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=PlaneStrain --generator=StandardSmallStrainMechanics --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=PlaneStress --generator=StandardFiniteStrainMechanics --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=PlaneStress --generator=StandardSmallStrainMechanics --source-file --header-file
#
if command -v clang-format &> /dev/null
then
  clang-format -i *xx
fi
# copying files
mv -- *BehaviourIntegrator.hxx ../include/MFEMMGIS/
mv -- *BehaviourIntegrator.cxx ../src

