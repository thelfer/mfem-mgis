/*!
 * \file   src/Material.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  Material::Material(std::shared_ptr<const PartialQuadratureSpace> s,
                     std::shared_ptr<const Behaviour> b_ptr)
      : MaterialDataManager(*b_ptr, s->getNumberOfIntegrationPoints()),
        quadrature_space(s),
        macroscopic_gradients(this->s1.gradients_stride, real(0)),
        behaviour_ptr(b_ptr) {}  // end of Material::Material

  void Material::setMacroscopicGradients(mgis::span<const real> g) {
    if (g.size() != this->s1.gradients_stride) {
      mgis::raise(
          "Material::setMacroscopicGradients: invalid number of components of "
          "the gradients");
    }
    std::copy(g.begin(), g.end(), this->macroscopic_gradients.begin());
  }  // end of Material::setMacroscopicGradients

  Material::~Material() = default;

};  // end of namespace mfem_mgis
