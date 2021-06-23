#!/bin/sh
#SBATCH -J mfem
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-core=1
#SBATCH -A b171
#SBATCH -t 0:55:00
#SBATCH -o period.%j.%a.out
#SBATCH -e period.%j.%a.err
#SBATCH --exclusive

# chargement des modules
module purge 
source /home/glatu/mod_mgis.sh

ulimit -a

MFEMMGIS_DIR=${HOME}/mfem-mgis
cd ${MFEMMGIS_DIR}/build/tests
MPIOPT="-report-bindings  --map-by core -bind-to core"
COMMONOPT="--mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh --library ./libBehaviourTest.so --test-case 0 --linearsolver 1 --refine 4"
time mpirun ./ssna303 --use-petsc --petsc-configuration-file rc_ex10p_jfnk
