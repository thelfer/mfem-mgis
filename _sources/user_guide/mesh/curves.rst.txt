===============
Defining curves
===============

This section describes the utility functions available in ``MFEM/MGIS``
for creating curves from parameters.

.. contents::
     :depth: 3
     :local:


Use Cases
^^^^^^^^^

The geometry utilities are particularly useful for:

- Defining boundary conditions or loads along specific paths in a mesh.
- Creating custom geometric features for post-processing or visualization.
- Generating points for interpolation or sampling along curves.

Overview
^^^^^^^^

The ``makePointsOnCurve`` function generates a set of points along a curve
based on a geometric description and a discretization strategy.

+--------------------------+-------------------------------------------------------------------+
| Function                 | Description                                                       |
+==========================+===================================================================+
| ``makePointsOnCurve<2>`` | Generates points along a curve in 2D space.                       |
+--------------------------+-------------------------------------------------------------------+
| ``makePointsOnCurve<3>`` | Generates points along a curve in 3D space.                       |
+--------------------------+-------------------------------------------------------------------+

**Parameters:**

The function accepts a ``Parameters`` object structured as follows:

- **Curve Type**: Currently, only the ``Line`` type is supported.

  - ``InitialPoint``: A vector of real numbers representing the starting point of the line.
  - ``FinalPoint``: A vector of real numbers representing the ending point of the line.
  - ``Discretization``: A nested ``Parameters`` object describing how the line is discretized.

**Discretization Strategies:**

1. **Uniform Discretization**

   - ``NumberOfPoints``: An integer specifying the number of points to generate along the line.
     Must be greater than or equal to 2.

   **Example:**

   .. code-block:: c++

      const auto params = mfem_mgis::Parameters{
          {"Line",
           mfem_mgis::Parameters{
               {"InitialPoint", std::vector<Parameter>{-2, 0.5}},
               {"FinalPoint", std::vector<Parameter>{6, 8}},
               {"Discretization",
                mfem_mgis::Parameters{
                    {"Uniform", mfem_mgis::Parameters{{"NumberOfPoints", 10}}}}}
           }}};
      const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, params);

2. **Geometric Progression Discretization**

   - ``NumberOfPoints``: An integer specifying the number of points to generate along the line.
     Must be greater than or equal to 2.
   - ``InitialDensity``: A positive real number representing the density of points at the start of the line.
   - ``FinalDensity``: A positive real number representing the density of points at the end of the line.

   **Example:**

   .. code-block:: c++

      const auto params = mfem_mgis::Parameters{
          {"Line",
           mfem_mgis::Parameters{
               {"InitialPoint", std::vector<Parameter>{2, 0}},
               {"FinalPoint", std::vector<Parameter>{17, 0}},
               {"Discretization",
                mfem_mgis::Parameters{
                    {"GeometricProgression",
                     mfem_mgis::Parameters{
                         {"NumberOfPoints", 10},
                         {"InitialDensity", real{1} / 150},
                         {"FinalDensity", real{1} / 3}}}}
           }}};
      const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, params);

.. note::

   The geometric progression discretization generates points with varying
   densities along the line, allowing for finer discretization at specific
   regions of the curve.

Error Handling
^^^^^^^^^^^^^^

Both ``makePoint`` and ``makePointsOnCurve`` functions return an
``std::optional`` type. If the input parameters are invalid or inconsistent
(e.g., wrong number of coordinates, missing parameters, invalid types),
the functions return an empty optional, and an error message is registered
in the context.

Common errors include:

- Missing required parameters (e.g., ``InitialPoint``, ``FinalPoint``, ``Discretization``).
- Incorrect parameter types (e.g., passing a string where a number is expected).
- Inconsistent space dimensions (e.g., providing 3 coordinates for a 2D point).
- Invalid discretization parameters (e.g., ``NumberOfPoints`` less than 2).

**Example of Error Handling:**

.. code-block:: c++

   const auto invalid_params = mfem_mgis::Parameters{
       {"Line",
        mfem_mgis::Parameters{
            {"InitialPoint", std::vector<Parameter>{-2, 0.5, 2}},  // Invalid: 3D point for 2D curve
            {"FinalPoint", std::vector<Parameter>{6, 8}},
            {"Discretization",
             mfem_mgis::Parameters{{"Uniform", mfem_mgis::Parameters{{"NumberOfPoints", 10}}}}
        }}};
   const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, invalid_params);
   if (!points) {
       // Handle error: points is empty
   }
