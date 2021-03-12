#!/bin/bash
cd ${HOME}/Work/Dev/mfem-mgis/build/tests
time mpirun -n 2 PeriodicTestP --mesh /home/latu/Work/Dev/mfem-mgis/tests/cube_2mat_per.mesh --library /home/latu/Work/Dev/mfem-mgis/build/tests/libBehaviourTest.so --test-case 0 --linearsolver 1
time mpirun -n 2 PeriodicTestNL --mesh /home/latu/Work/Dev/mfem-mgis/tests/cube_2mat_per.mesh --library /home/latu/Work/Dev/mfem-mgis/build/tests/libBehaviourTest.so --test-case 0 --algo 0 --linearsolver 1
time mpirun -n 2 PeriodicTestNL --mesh /home/latu/Work/Dev/mfem-mgis/tests/cube_2mat_per.mesh --library /home/latu/Work/Dev/mfem-mgis/build/tests/libBehaviourTest.so --test-case 0 --algo 1 --linearsolver 1
