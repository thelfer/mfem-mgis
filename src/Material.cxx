/*!
 * \file   src/Material.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  Material::Material(std::shared_ptr<const PartialQuadratureSpace> s,
                     std::shared_ptr<const Behaviour> b_ptr)
      : MaterialDataManager(*b_ptr, s->getNumberOfIntegrationPoints()),
        quadrature_space(s),
        behaviour_ptr(b_ptr) {}  // end of Material::Material

  Material::~Material() = default;

};  // end of namespace mfem_mgis

