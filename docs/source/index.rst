.. _mfem_mgis:

=========================
The ``MFEM/MGIS`` project
=========================

.. toctree::
   :hidden:
   :maxdepth: 1

   user_guide/user_guide.rst
   commented_examples/commented_examples.rst
   developer_guide/developer_guide.rst
   installation_guide/installation_guide.rst

The aim of the ``MFEM/MGIS`` project is to provide a ``C++`` library to
build advanced mechanical simulation based on ``MFront`` behaviours with
``MFEM`` . The ``MGIS`` project is used to interface
``MFront`` behaviours .

The ``MFEM-MGIS`` project, aims at efficiently use supercomputers in
order to describe coupled multiphysics phenomena with a particular focus
on thermo-mechanics. This open-source library is based on several
components as prerequisites:

- the ``mfem`` (Modular Finite Element Methods) library
  :cite:`Anderson_2020`,
- the ``MGIS`` (MFront Generic Interface Support) library
  :cite:`helfer_mfrontgenericinterfacesupport_2020`,
-  the ``MFront`` code generator.

Thanks to the features embedded within ``mgis`` and ``MFront`` and
thanks to specific developments, ``MFEM-MGIS`` adds several mechanical
features compared to a pure ``mfem`` approach.

The library tackles some peculiarities of nonlinear mechanics. In
particular, the support of complex constitutive laws and the management
of advanced boundary conditions. It provides a high level of abstraction
based focused on the physics to be treated.

.. figure:: img/mfem-mgis-illustration.png

References
==========

.. bibliography:: bibliography.bib
