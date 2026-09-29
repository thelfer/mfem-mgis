#!/bin/bash

sed \
    -e "s%mesh_file = nullptr%mesh_file = \"cube.mesh\"%g" \
    -e "s%behaviour = nullptr%behaviour = \"IsotropicLinearHardeningPlasticity\"%g" \
    -e "s%library = nullptr%library = \"src/libBehaviour.so\"%g" \
    -e "s%isv_name = nullptr%isv_name = \"EquivalentPlasticStrain\"%g" \
    "$1" > "$2"
    
