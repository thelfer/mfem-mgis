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
# non linear heat transfer
./behaviour-integrator -s isotropic --hypothesis=PlaneStress --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file
./behaviour-integrator -s orthotropic --hypothesis=PlaneStress --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=PlaneStrain --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file
./behaviour-integrator -s orthotropic --hypothesis=PlaneStrain --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file
./behaviour-integrator -s isotropic --hypothesis=Tridimensional --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file
./behaviour-integrator -s orthotropic --hypothesis=Tridimensional --generator=StationaryNonLinearHeatTransfer --unknown-name=Temperature --requires-unknowns-values --source-file --header-file

#
if command -v clang-format &> /dev/null
then
  clang-format -i *xx
fi
# copying files
mv -- *BehaviourIntegrator.hxx ../include/MFEMMGIS/
mv -- *BehaviourIntegrator.cxx ../src

