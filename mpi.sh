#!/bin/bash
MFEMMGIS_DIR=${HOME}/mfem-mgis
cd ${MFEMMGIS_DIR}/build/tests
time mpirun -n 4 PeriodicTestP --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh --library ${MFEMMGIS_DIR}/build/tests/libBehaviourTest.so --test-case 0 --linearsolver 1
time mpirun -n 4 PeriodicTestNL --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh --library ${MFEMMGIS_DIR}/build/tests/libBehaviourTest.so --test-case 0 --algo 0 --linearsolver 1
time mpirun -n 4 PeriodicTestNL --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh --library ${MFEMMGIS_DIR}/build/tests/libBehaviourTest.so --test-case 0 --algo 1 --linearsolver 1
