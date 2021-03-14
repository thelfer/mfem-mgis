MGIS_PACKAGE=/jjivefi5pmeimwn2
cmake .. -DCMAKE_BUILD_TYPE=Release  \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_INSTALL_PREFIX=$PWD/../install \
    -DMFrontGenericInterface_DIR=$(spack location -i ${MGIS_PACKAGE})/share/mgis/cmake
make -j
