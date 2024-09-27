/*!
 * \file   include/MFEMMGIS/TridimensionalUniformImposedPressureBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_TRIDIMENSIONALUNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_TRIDIMENSIONALUNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX

#include <memory>
#include <mfem/linalg/densemat.hpp>
#include "mfem/fem/nonlininteg.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT TridimensionalUniformImposedPressureBoundaryCondition final
      : public NonlinearFormIntegrator {
    /*!
     * \brief constructor
     * \param[in] prvalues: function returning the imposed pressure
     */
    TridimensionalUniformImposedPressureBoundaryCondition(
        std::function<real(const real)>);
    // MFEM API
    void AssembleElementVector(const mfem::FiniteElement &,
                               mfem::ElementTransformation &,
                               const mfem::Vector &,
                               mfem::Vector &) override;

    void AssembleElementGrad(const mfem::FiniteElement &,
                             mfem::ElementTransformation &,
                             const mfem::Vector &,
                             mfem::DenseMatrix &) override;
    //! \brief destructor
    virtual ~TridimensionalUniformImposedPressureBoundaryCondition();

   protected:
    //! \brief function returning the value of the imposed pressure
    std::function<real(const real)> prfct;

#ifndef MFEM_THREAD_SAFE
   private:
    //! \brief vector used to store the value of the shape functions
    mfem::Vector shape;
#endif
  };  // end of TridimensionalUniformImposedPressureBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_TRIDIMENSIONALUNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX */
