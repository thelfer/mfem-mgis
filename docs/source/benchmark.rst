Benchmarks
==========


Benchmarks REV MOX
^^^^^^^^^^^^^^^^^^

In this section, we propose some scaling curves for different test cases. For the first benchmark, we use the MOX REV example with 643 inclusions and a viscoplastic behavior law for the matrix, and an elastic behavior law for the inclusions. Calculations are performed on the TOPAZE supercalculator at CCRT. Each node in the cluster is built on 64-core AMD EPYC Milan 7763 dual-socket processors running at 2.45 GHz and equipped with 256 GB RAM. Benchmarks are in pure MPI. 


Regarding the specificities of the simulations, we use the HyprePCG solver with a HypreBoomerAMG preconditioner. Mesh reading is performed using a mesh pre-cut into small msh files to limit the impact on the memory footprint.


.. figure:: _static/634.jpeg
   :alt: Illustration of a RVE with 634 spheres after 5 seconds.


We performed tests on 3 problem sizes: 80M ddl, 190M ddl and 664M ddl.


.. note::

  Un test de faisabilité de cette simulation a aussi était réalisé sûr plus de 5.3 milliards de ddl.


80M ddl
-------

.. figure:: _static/80MDofMFEMMGIS.png
   :alt: Time, Memory footprint and speedup of a MOX RVE with 80M ddl.


190M ddl
--------

.. figure:: _static/190MDofMFEM-MGIS.png
   :alt: Time, Memory footprint and speedup of a MOX RVE with 190M ddl.

664M ddl
---------

.. figure:: _static/664MDofMFEM-MGIS.png
   :alt: Time, Memory footprint and speedup of a MOX REV with 664M ddl.


MFEM-MGIS / MFEM 
^^^^^^^^^^^^^^^^

This benchmark aims to evaluate the time overhead between MFEM and MFEM-MGIS for a simple test case involving an elastic behavior law. The same supercomputer used in the previous RVE benchmark is employed.

.. figure:: _static/bench-perf.png



+------------+-----+---------+----------------+---------------------+------------+----------------+----------+
| Refinement | MPI |   DDL   | MFEM-MGIS Time | MFEM-MGIS Iteration | MFEM Time  | MFEM Iteration | Overhead |
+------------+-----+---------+----------------+---------------------+------------+----------------+----------+
| 4          |  1  | 111843  |          8.34  |                1592 |       6.21 |           1565 |  34.4 %  |
+------------+-----+---------+----------------+---------------------+------------+----------------+----------+
| 5          |  1  | 839619  |        122.54  |                3167 |      94.53 |           3102 |  29.6 %  |
+------------+-----+---------+----------------+---------------------+------------+----------------+----------+
| 5          | 16  | 839619  |         51.34  |                3167 |      41.90 |           3102 |  22.5 %  |
+------------+-----+---------+----------------+---------------------+------------+----------------+----------+
| 6          | 16  | 6502275 |       1057.24  |                6273 |     877.92 |           6124 |  20.0 %  |
+------------+-----+---------+----------------+---------------------+------------+----------------+----------+


.. note::

  Note that due to the small size of the input mesh, domain decomposition is limited to 16 subdomains (MFEM-MGIS only performs parallel refinement). To improve this benchmark, consider using sequential refinement or a finer mesh.
  
Run this example on Topaze supercomputer:

.. code-block:: bash

  ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5
  ccc_mprun -n 16 -c 1 -m work -T 600 -p milan -Q test ./MFEMMGISLinearElasticityBenchmark --mesh ../beam-tet.mesh -r 5 

