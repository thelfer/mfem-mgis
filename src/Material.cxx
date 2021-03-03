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

  [[noreturn]] static std::array<real, 9u> raiseInvalidGetRotationMatrixCall(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    mgis::raise(
        "Material::getRotationMatrix: "
        "invalid call for isotropic behaviours");
  }  // end of raiseInvalidGetRotationMatrixCall

  [[noreturn]] static std::array<real, 9u> raiseUnsetRotationMatrix(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    mgis::raise(
        "Material::getRotationMatrix: "
        "unset rotation matrix");
  }  // end of raiseInvalidGetRotationMatrixCall

  Material::Material(std::shared_ptr<const PartialQuadratureSpace> s,
                     std::unique_ptr<const Behaviour> b_ptr)
      : MaterialDataManager(*b_ptr, s->getNumberOfIntegrationPoints()),
        quadrature_space(s),
        macroscopic_gradients(this->s1.gradients_stride, real(0)),
        behaviour_ptr(std::move(b_ptr)),
        get_rotation_fct_ptr(
            this->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC
                ? raiseInvalidGetRotationMatrixCall
                : raiseUnsetRotationMatrix) {}  // end of Material::Material

  void Material::setMacroscopicGradients(mgis::span<const real> g) {
    if (g.size() != this->s1.gradients_stride) {
      mgis::raise(
          "Material::setMacroscopicGradients: "
          "invalid number of components of the gradients");
    }
    std::copy(g.begin(), g.end(), this->macroscopic_gradients.begin());
  }  // end of setMacroscopicGradients

  static void checkBehaviourSymmetry(const Behaviour &b) {
    if (b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      mgis::raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour is not orthotropic)");
    }
  }  // end of checkBehaviourSymmetry

  void Material::setRotationMatrix(const RotationMatrix2D &r) {
    checkBehaviourSymmetry(b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 2u) {
      mgis::raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour' hypothesis is not 2D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &rm,  //
              const RotationMatrix3D &,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else {
      mgis::raise("Material::getRotationMatrix: unimplemented case yet");
    }
    this->r2D = r;
  }  // end of setRotationMatrix

  void Material::setRotationMatrix(const RotationMatrix3D &r) {
    checkBehaviourSymmetry(b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 3u) {
      mgis::raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour hypothesis is not 3D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &,  //
              const RotationMatrix3D &rm,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else {
      mgis::raise("Material::getRotationMatrix: unimplemented case yet");
    }
    this->r3D = r;
  }  // end of setRotationMatrix

  //   std::array<real, 9u> Material::getRotationMatrix(const size_type) const {
  //     if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) == 3u) {
  //       if (std::holds_alternative<std::array<real, 9u>>(this->r3D)) {
  //         return std::get<std::array<real, 9u>>(this->r3D);
  //       } else {
  //         const auto &axes = std::get<std::array<MaterialAxis3D,
  //         2u>>(this->r3D); if ((std::holds_alternative<std::array<real,
  //         3u>>(axes[0])) &&
  //             (std::holds_alternative<std::array<real, 3u>>(axes[1]))) {
  //           const auto &a1 = std::get<std::array<real, 3u>>(axes[0]);
  //           const auto &a2 = std::get<std::array<real, 3u>>(axes[1]);
  //           return mgis::behaviour::buildRotationMatrix(a1, a2);
  //         } else {
  //           mgis::raise("Material::getRotationMatrix: unimplemented case
  //           yet");
  //         }
  //       }
  //     } else if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) ==
  //     2u) {
  //       if (std::holds_alternative<std::array<real, 9u>>(this->r2D)) {
  //         return std::get<std::array<real, 9u>>(this->r2D);
  //       } else {
  //         mgis::raise("Material::getRotationMatrix: unimplemented case yet");
  //       }
  //     } else if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) !=
  //     1u) {
  //       mgis::raise("Material::getRotationMatrix: unimplemented case yet");
  //     }
  //     return {1, 0, 0,  //
  //             0, 1, 0,  //
  //             0, 0, 1};
  //   }  // end of Material::getRotationMatrix

  Material::~Material() = default;

};  // end of namespace mfem_mgis
