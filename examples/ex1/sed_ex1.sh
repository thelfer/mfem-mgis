#!/bin/bash

sed \
    -e "s%mesh_file = nullptr%mesh_file = \"cube.mesh\"%g" \
    -e "s%behaviour = nullptr%behaviour = \"Plasticity\"%g" \
    -e "s%library = nullptr%library = \"src/libBehaviour.so\"%g" \
    -e "s%reference_file = nullptr%reference_file = \"Plasticity.ref\"%g" \
    -e "s%isv_name = nullptr%isv_name = \"EquivalentPlasticStrain\"%g" \
    "$1" > "$2"
    
