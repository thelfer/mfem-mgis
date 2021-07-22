/*!
 * \file   src/Material.cxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/RotationMatrix.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  [[noreturn]] static std::array<real, 9u> raiseInvalidGetRotationMatrixCall(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    raise(
        "Material::getRotationMatrix: "
        "invalid call for isotropic behaviours");
  }  // end of raiseInvalidGetRotationMatrixCall

  [[noreturn]] static std::array<real, 9u> raiseUnsetRotationMatrix(
      const RotationMatrix2D &, const RotationMatrix3D &, const size_type) {
    raise(
        "Material::getRotationMatrix: "
        "unset rotation matrix");
  }  // end of raiseInvalidGetRotationMatrixCall

  Material::Material(std::shared_ptr<const PartialQuadratureSpace> s,
                     std::unique_ptr<const Behaviour> b_ptr)
      : MaterialDataManager(*b_ptr, s->getNumberOfIntegrationPoints()),
        quadrature_space(s),
        macroscopic_gradients(this->s1.gradients_stride, real(0)),
        get_rotation_fct_ptr(this->b.symmetry ==
                                     mgis::behaviour::Behaviour::ORTHOTROPIC
                                 ? &raiseInvalidGetRotationMatrixCall
                                 : &raiseUnsetRotationMatrix),
        behaviour_ptr(std::move(b_ptr)) {
    this->allocateArrayOfTangentOperatorBlocks();
}  // end of Material::Material

  const PartialQuadratureSpace &Material::getPartialQuadratureSpace() const {
    return *(this->quadrature_space);
  }  // end of getPartialQuadratureSpace

  std::shared_ptr<const PartialQuadratureSpace>
  Material::getPartialQuadratureSpacePointer() const {
    return this->quadrature_space;
  }  // end of getPartialQuadratureSpacePointer

  void Material::setMacroscopicGradients(mgis::span<const real> g) {
    if (g.size() != this->s1.gradients_stride) {
      raise(
          "Material::setMacroscopicGradients: "
          "invalid number of components of the gradients");
    }
    std::copy(g.begin(), g.end(), this->macroscopic_gradients.begin());
  }  // end of setMacroscopicGradients

  static void checkBehaviourSymmetry(const Behaviour &b) {
    if (b.symmetry != mgis::behaviour::Behaviour::ORTHOTROPIC) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour is not orthotropic)");
    }
  }  // end of checkBehaviourSymmetry

  void Material::setRotationMatrix(const RotationMatrix2D &r) {
    checkBehaviourSymmetry(b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 2u) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour' hypothesis is not 2D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &rm, const RotationMatrix3D &,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else if (std::holds_alternative<
                   std::shared_ptr<PartialQuadratureFunction>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &rm, const RotationMatrix3D &,
              const size_type o) {
            std::array<real, 9u> m;
            const auto &f =
                *(std::get<std::shared_ptr<PartialQuadratureFunction>>(rm));
            return mgis::behaviour::buildRotationMatrix(
                f.getIntegrationPointValues<2>(o));
            return m;
          };
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r2D = r;
  }  // end of setRotationMatrix

  void Material::setRotationMatrix(const RotationMatrix3D &r) {
    checkBehaviourSymmetry(b);
    if (mgis::behaviour::getSpaceDimension(this->b.hypothesis) != 3u) {
      raise(
          "Material::setRotationMatrix: "
          "invalid call (behaviour hypothesis is not 3D)");
    }
    if (std::holds_alternative<std::array<real, 9u>>(r)) {
      this->get_rotation_fct_ptr =
          +[](const RotationMatrix2D &,  //
              const RotationMatrix3D &rm,
              const size_type) { return std::get<std::array<real, 9u>>(rm); };
    } else if (std::holds_alternative<std::array<MaterialAxis3D, 2u>>(r)) {
      const auto &[v1, v2] = std::get<std::array<MaterialAxis3D, 2u>>(r);
      if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
          (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
                 (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<std::array<real, 3u>>(v1)) &&
                 (std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f2 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a2));
          return mgis::behaviour::buildRotationMatrix(
              std::get<std::array<real, 3u>>(a1),
              f2.getIntegrationPointValues<3>(o));
        };
      } else if ((std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v1)) &&
                 (std::holds_alternative<std::array<real, 3u>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f1 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a1));
          return mgis::behaviour::buildRotationMatrix(
              f1.getIntegrationPointValues<3>(o),
              std::get<std::array<real, 3u>>(a2));
        };
      } else if ((std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v1)) &&
                 (std::holds_alternative<
                     std::shared_ptr<PartialQuadratureFunction>>(v2))) {
        this->get_rotation_fct_ptr = +[](const RotationMatrix2D &,  //
                                         const RotationMatrix3D &rm,
                                         const size_type o) {
          const auto &[a1, a2] = std::get<std::array<MaterialAxis3D, 2u>>(rm);
          const auto &f1 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a1));
          const auto &f2 =
              *(std::get<std::shared_ptr<PartialQuadratureFunction>>(a2));
          return mgis::behaviour::buildRotationMatrix(
              f1.getIntegrationPointValues<3>(o),
              f2.getIntegrationPointValues<3>(o));
        };
      } else {
        raise("Material::setRotationMatrix: unimplemented case yet");
      }
    } else {
      raise("Material::setRotationMatrix: unimplemented case yet");
    }
    this->r3D = r;
  }  // end of setRotationMatrix

  Material::~Material() = default;

};  // end of namespace mfem_mgis
