#! /usr/bin/env bash
set -e

# compilation of the code generator
g++ -Wall -W -pedantic -Werror -std=c++17 -g behaviour-integrator.cxx -o behaviour-integrator -lginac -lcln
#
for hypothesis in Tridimensional PlaneStrain PlaneStress;
do
  for symmetry in isotropic orthotropic
  do
    for generator in StandardFiniteStrainMechanics StandardSmallStrainMechanics
    do
      ./behaviour-integrator -s ${symmetry} --hypothesis=${hypothesis} --generator=${generator} --source-file --header-file
    done
  done
done
#
if command -v clang-format &> /dev/null
then
  clang-format -i *xx
fi
# copying files
mv -- *BehaviourIntegrator.hxx ../include/MFEMMGIS/
mv -- *BehaviourIntegrator.cxx ../src

