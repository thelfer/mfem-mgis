#!/bin/bash

sed "s%mesh_file = nullptr%mesh_file = \"cube.mesh\"%g;s%behaviour = nullptr%behaviour = \"Plasticity\"%g;s%library = nullptr%library = \"src/libBehaviour.so\"%g;s%reference_file = nullptr%reference_file = \"Plasticity.ref\"%g;s%isv_name = nullptr%isv_name = \"EquivalentPlasticStrain\"%g" $1 > $2

