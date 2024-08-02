=============================================================
Implementing a micromorphic damage behaviour in ``MFEM/MGIS``
=============================================================

:Author: Thomas Helfer
:Date:   2021

Implementation of a dedicated behaviour integrator
==================================================

Basic header file
-----------------

.. code:: cpp

   /*!
    * \file   include/MFEMMGIS/MicromorphicDamage2DBehaviourIntegrator.hxx
    * \brief
    * \author Thomas Helfer
    * \date   07/12/2021
    * \brief header file declaring the MicromorphicDamage2DBehaviourIntegrator
    * class
    */

   #ifndef LIB_MFEM_MGIS_MICROMORPHICDAMAGE2DBEHAVIOURINTEGRATOR_HXX
   #define LIB_MFEM_MGIS_MICROMORPHICDAMAGE2DBEHAVIOURINTEGRATOR_HXX

   #include "MFEMMGIS/BehaviourIntegratorBase.hxx"

   namespace mfem_mgis {

     // forward declaration
     struct FiniteElementDiscretization;

     /*!
      * \brief class implementing the a behaviour integrator dedicated to
      * micromorphic damage in two dimensions (the modelling hypothesis has no
      * effect on this specific behaviour as out of plane damage gradients are
      * assumed to be zero).
      */
     struct MFEM_MGIS_EXPORT MicromorphicDamage2DBehaviourIntegrator
         : BehaviourIntegratorBase {
       /*!
        * \brief constructor
        * \param[in] s: quadrature space
        * \param[in] b_ptr: behaviour
        */
       MicromorphicDamage2DBehaviourIntegrator(
           std::shared_ptr<const PartialQuadratureSpace>,
           std::unique_ptr<const Behaviour>);
       /*!
        * \brief constructor
        * \param[in] s: quadrature space
        * \param[in] m: material attribute.
        * \param[in] b_ptr: behaviour
        */
       MicromorphicDamage2DBehaviourIntegrator(
           const FiniteElementDiscretization &,
           const size_type,
           std::unique_ptr<const Behaviour>);
       //
       const mfem::IntegrationRule &getIntegrationRule(
           const mfem::FiniteElement &,
           const mfem::ElementTransformation &) const override;
       real getIntegrationPointWeight(mfem::ElementTransformation &,
                                      const mfem::IntegrationPoint &) const
           noexcept override;
       bool integrate(const mfem::FiniteElement &,
                      mfem::ElementTransformation &,
                      const mfem::Vector &,
                      const IntegrationType) override;

       void updateResidual(mfem::Vector &,
                           const mfem::FiniteElement &,
                           mfem::ElementTransformation &,
                           const mfem::Vector &) override;

       void updateJacobian(mfem::DenseMatrix &,
                           const mfem::FiniteElement &,
                           mfem::ElementTransformation &,
                           const mfem::Vector &) override;

       void computeInnerForces(mfem::Vector &,
                               const mfem::FiniteElement &,
                               mfem::ElementTransformation &) override;
       //! \brief destructor
       ~MicromorphicDamage2DBehaviourIntegrator() override;

      private:
       /*!
        * \return the integration rule for the given element and  * element
        * transformation. \param[in] e: element \param[in] tr: element
        * transformation
        */
       static const mfem::IntegrationRule &selectIntegrationRule(
           const mfem::FiniteElement &, const mfem::ElementTransformation &);
       /*!
        * \brief build the quadrature space for the given  * material
        * \param[in] fed: finite element discretization.
        * \param[in] m: material attribute.
        */
       static std::shared_ptr<const PartialQuadratureSpace> buildQuadratureSpace(
           const FiniteElementDiscretization &, const size_type);

   #ifndef MFEM_THREAD_SAFE
       //! \brief vector used to store the value of the shape functions
       mfem::Vector shape;
       //! \brief matrix used to store the derivatives of the shape functions
       mfem::DenseMatrix dshape;
   #endif
     }; // end of struct MicromorphicDamage2DBehaviourIntegrator

   } // end of namespace mfem_mgis

   #endif /* LIB_MFEM_MGIS_MICROMORPHICDAMAGE2DBEHAVIOURINTEGRATOR_HXX */

Implementation
--------------

Implementation of the ``selectIntegrationRule`` static method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

     const mfem::IntegrationRule &
     MicromorphicDamage2DBehaviourIntegrator::selectIntegrationRule(
         const mfem::FiniteElement &e, const mfem::ElementTransformation &t) {
       const auto order = 2 * t.OrderGrad(&e);
       return mfem::IntRules.Get(e.GetGeomType(), order);
     }

Implementation of the ``buildQuadratureSpace`` static method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

     std::shared_ptr<const PartialQuadratureSpace>
     MicromorphicDamage2DBehaviourIntegrator::buildQuadratureSpace(
         const FiniteElementDiscretization &fed, const size_type m) {
       auto selector = [](const mfem::FiniteElement &e,
                          const mfem::ElementTransformation &tr)
           -> const mfem::IntegrationRule & {
         return selectIntegrationRule(e, tr);
       };  // end of selector
       return std::make_shared<PartialQuadratureSpace>(fed, m, selector);
     }  // end of buildQuadratureSpace

Implementation of ``MicromorphicDamage2DBehaviourIntegrator`` constructor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

     MicromorphicDamage2DBehaviourIntegrator::
         MicromorphicDamage2DBehaviourIntegrator(
             const FiniteElementDiscretization &fed,
             const size_type m,
             std::unique_ptr<const Behaviour> b_ptr)
         : BehaviourIntegratorBase(buildQuadratureSpace(fed, m),
                                   std::move(b_ptr)) {
       if (this->b.symmetry != Behaviour::ISOTROPIC) {
         raise("invalid behaviour symmetry");
       }
     }  // end of MicromorphicDamage2DBehaviourIntegrator

Implementation of the ``getIntegrationPointWeight`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

     real MicromorphicDamage2DBehaviourIntegrator::getIntegrationPointWeight(
         mfem::ElementTransformation &tr, const mfem::IntegrationPoint &ip) const
         noexcept {
       return ip.weight * tr.Weight();
     } // end of getIntegrationPointWeight

Implementation of the ``getIntegrationRule`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: cpp

     const mfem::IntegrationRule &
     MicromorphicDamage2DBehaviourIntegrator::getIntegrationRule(
         const mfem::FiniteElement &e,
         const mfem::ElementTransformation &t) const {
       return MicromorphicDamage2DBehaviourIntegrator::selectIntegrationRule(e, t);
     }  // end of getIntegrationRule

Modifying ``CMakeLists.txt`` files
----------------------------------

Declaration in the ``BehaviourIntegrator`` factory
--------------------------------------------------

.. code:: cpp

   #include "MFEMMGIS/MicromorphicDamage2DBehaviourIntegrator.hxx"

.. code:: cpp

       f.addGenerator(
           "MicromorphicDamage",  //
           [](const FiniteElementDiscretization& fed, const size_type m,
              std::unique_ptr<const Behaviour> b) {
             return std::make_unique<MicromorphicDamage2DBehaviourIntegrator>(
                 fed, m, std::move(b));
           });
