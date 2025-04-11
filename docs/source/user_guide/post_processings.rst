.. _mfem_mgis_post_processings:

================
Post-processings
================

In this section, we describe the post-processing options available in
MFEM-MGIS.

.. contents::
    :depth: 3
    :local:

.. warning::

  This section is under construction

The :code:`ParaviewExport` post-processing
==========================================

This post-processing allows to export the unknowns of a nonlinear
evolution problem for visualization in :code:`paraview`:

- Key: ``ParaviewExportResults``

**Example:**

.. code-block:: cpp

  problem.addPostProcessing("ParaviewExportResults",
                            {{"OutputFileName", "SatohTestOutput"}});

**Results**

.. figure:: img/SatohTest.png

It is also possible to extract only portions of the mesh by defining either the boundary zones or the domain attribute. This is particularly useful for reducing the size of output files when only a small part is to be studied.

**Example**

.. code-block:: cpp

    std::vector<int> DomainAttibuteLeft = {1};
    std::vector<int> DomainAttibuteRight = {2};
    std::vector<int> AllBoundaries = {1, 2};
    /** You can not define DomainAttributes and BoundaryAttributes in a single post processing */
    problem.addPostProcessing("ParaviewExportResults",
        {{"OutputFileName", "TestPPSubMeshOutputDir/AllMesh"},
        {"OutputFieldName", "Displacement"},
        {"Verbosity", 1}});
    problem.addPostProcessing("ParaviewExportResults",
        {{"OutputFileName", "TestPPSubMeshOutputDir/Attribute1"},
        {"OutputFieldName", "Displacement"},
        {"DomainAttributes", DomainAttibuteLeft},
        {"Verbosity", 1}});
    problem.addPostProcessing("ParaviewExportResults",
        {{"OutputFileName", "TestPPSubMeshOutputDir/Attribute2"},
        {"OutputFieldName", "Displacement"},
        {"DomainAttributes", DomainAttibuteRight},
        {"Verbosity", 1}});
    problem.addPostProcessing("ParaviewExportResults",
        {{"OutputFileName", "TestPPSubMeshOutputDir/Boundaries"},
        {"OutputFieldName", "Displacement"},
        {"BoundaryAttributes", AllBoundaries},
        {"Verbosity", 1}});

**Results**

.. figure:: img/ExportDataAtNodes.png


+---------------------+--------------------------------------------------------------------------------------------+
| **Key**             | **Description**                                                                            |
+=====================+============================================================================================+
| OutputFileName      | Name of the output directory                                                               |
+---------------------+--------------------------------------------------------------------------------------------+
| OutputFieldName     | Name of the field that will appear in ParaView (default: ``"u"``)                          |
+---------------------+--------------------------------------------------------------------------------------------+
| DomainAttributes    | List of domain attributes; a submesh will be used instead of exporting the entire mesh     |
+---------------------+--------------------------------------------------------------------------------------------+
| BoundaryAttributes  | List of boundary attributes; a submesh will be used instead of exporting the entire mesh   |
+---------------------+--------------------------------------------------------------------------------------------+
| Verbosity           | If this value is ``>= 1``, submesh information will be displayed when using attributes     |
+---------------------+--------------------------------------------------------------------------------------------+

Export Integration Point Results At Nodes
==========================================

- Key: ``ParaviewExportIntegrationPointResultsAtNodes``

**Example:**

.. code-block:: cpp

  auto results = std::vector<mfem_mgis::Parameter>{
      "Stress", "ImposedTemperature", "HydrostaticPressure"};
  problem.addPostProcessing(
      "ParaviewExportIntegrationPointResultsAtNodes",
      {{"OutputFileName", "SatohTestIntegrationPointOutput"},
       {"Materials", {"plate"}},
       {"Results", results}});

**Results**

.. figure:: img/SatohTestStress.png

.. figure:: img/SatohTestTemperature.png

.. figure:: img/SatohTestPressure.png


Compute Mean Thermodynamic Forces
=================================

 The `Compute Mean Thermodynamic Forces` post-processing step calculates the average stress over selected regions of the mesh. 

- Key: ``MeanThermodynamicForces``

**Example: print the average stress of one inclusion into a matrix (RVE)**

.. code-block:: cpp

  p.addPostProcessing(
      "MeanThermodynamicForces",
      {{"OutputFileName", "avgStress"}});

**Results**

We display the average stress SZZ over the RVE (composed of 83% matrix and 17% inclusion), we process the avgstress file and then plot the result: 

.. code-block:: text

  awk '{if(NR>13) print $1 " " 0.83*$4+0.17*$10}' avgStress > res-mfem-mgis.txt
  plot "res-mfem-mgis.txt" u 1:2 w l title "mfem-mgis"

.. figure:: img/avgStress.png


Compute Stored Energy
=====================

- Key: ``StoredEnergy``

**Example:**

.. code-block:: cpp

  p.addPostProcessing(
      "StoredEnergy",
      {{"OutputFileName", "energy.txt"}});

Compute dissipated Energy
=========================

- Key: ``DissipatedEnergy``

**Example:**

.. code-block:: cpp

  p.addPostProcessing(
      "DissipatedEnergy",
      {{"OutputFileName", "dissiped_energy.txt"}});

