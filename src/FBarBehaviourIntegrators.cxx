/*!
 * \file src/FBarBehaviourIntegrators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   22/03/2026
 */

#include "MFEMMGIS/FBarBehaviourIntegrators.hxx"
#include "MFEMMGIS/FBarIsotropicPlaneStrainBehaviourIntegrator.hxx"
#include "MFEMMGIS/FBarIsotropicTridimensionalBehaviourIntegrator.hxx"
#include "MFEMMGIS/FBarOrthotropicPlaneStrainBehaviourIntegrator.hxx"
#include "MFEMMGIS/FBarOrthotropicTridimensionalBehaviourIntegrator.hxx"

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
      return ctx.registerErrorMessage(
          "invalid behaviour type: a finite strain behaviour is expected");
    }
    if (b->symmetry == Behaviour::ISOTROPIC) {
      auto bi = make_unique<FBarIsotropicPlaneStrainBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
      if (isInvalid(bi)) {
        return {};
      }
      const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
      // getMaterial can't fail here
      bi->getMaterial(ctx)->setMacroscopicGradients(F);
      return bi;
    }
    auto bi = make_unique<FBarOrthotropicPlaneStrainBehaviourIntegrator>(
        ctx, fed, m, std::move(b));
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
    // getMaterial can't fail here
    bi->getMaterial(ctx)->setMacroscopicGradients(F);
    return bi;
  }  // end of generatePlaneStrainFBarBehaviourIntegrators

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
      return ctx.registerErrorMessage(
          "invalid behaviour type: a finite strain behaviour is expected");
    }
    if (b->symmetry == Behaviour::ISOTROPIC) {
      auto bi = make_unique<FBarIsotropicTridimensionalBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
      if (isInvalid(bi)) {
        return {};
      }
      const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
      // getMaterial can't fail here
      bi->getMaterial(ctx)->setMacroscopicGradients(F);
      return bi;
    }
    auto bi = make_unique<FBarOrthotropicTridimensionalBehaviourIntegrator>(
        ctx, fed, m, std::move(b));
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    // getMaterial can't fail here
    bi->getMaterial(ctx)->setMacroscopicGradients(F);
    return bi;
  }  // end of generateTridimensionalFBarBehaviourIntegrators

}  // end of namespace mfem_mgis
