/*!
 * \file src/FBarBehaviourIntegrators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   22/03/2026
 */

#include "MFEMMGIS/FBarBehaviourIntegrators.hxx"
#include "MFEMMGIS/FBarIsotropicPlaneStrainBehaviourIntegrator.hxx"
#include "MFEMMGIS/FBarIsotropicPlaneStressBehaviourIntegrator.hxx"
#include "MFEMMGIS/FBarIsotropicTridimensionalBehaviourIntegrator.hxx"

namespace mfem_mgis {

  std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFBarBehaviourIntegrators(
      Context &ctx,
      const FiniteElementDiscretization &fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters &params) noexcept {
    if (!params.empty()) {
      return ctx.registerErrorMessage("no parameter expected");
    }
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<FBarIsotropicPlaneStrainBehaviourIntegrator>(
        ctx, fed, m, std::move(b));
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generatePlaneStrainFBarBehaviourIntegrators

  std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressFBarBehaviourIntegrators(
      Context &ctx,
      const FiniteElementDiscretization &fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters &params) noexcept {
    if (!params.empty()) {
      return ctx.registerErrorMessage("no parameter expected");
    }
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<FBarIsotropicPlaneStressBehaviourIntegrator>(
        ctx, fed, m, std::move(b));
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generatePlaneStressFBarBehaviourIntegrators

  std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFBarBehaviourIntegrators(
      Context &ctx,
      const FiniteElementDiscretization &fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters &params) noexcept {
    if (!params.empty()) {
      return ctx.registerErrorMessage("no parameter expected");
    }
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<FBarIsotropicTridimensionalBehaviourIntegrator>(
        ctx, fed, m, std::move(b));
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generateTridimensionalFBarBehaviourIntegrators

}  // end of namespace mfem_mgis
