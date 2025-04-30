/*!
 * \file   src/PartialQuadratureFunctionEvaluator.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#include "MFEMMGIS/PartialQuadratureFunctionEvaluator.hxx"

namespace mfem_mgis{

  RotationMatrixPartialQuadratureFunctionEvalutor::
      RotationMatrixPartialQuadratureFunctionEvalutor(const Material& m)
      : material(m) {
    raise_if(
        this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC,
        "material is not orthotropic");
  }

  void RotationMatrixPartialQuadratureFunctionEvalutor::check() const {
    raise_if(
        this->material.b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC,
        "material is not orthotropic");
  } // end of check

  const PartialQuadratureSpace&
  RotationMatrixPartialQuadratureFunctionEvalutor::getPartialQuadratureSpace()
      const {

    return this->material.getPartialQuadratureSpace();
  }

} // end of namespace mfem_mgis
