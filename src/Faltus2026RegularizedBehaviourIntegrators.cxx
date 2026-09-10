/*!
 * \file   Faltus2026RegularizedBehaviourIntegrators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   17/03/2026
 */

#include <array>
#include "mfem/fem/fe.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.hxx"

namespace mfem_mgis {

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<Faltus2026RegularizedIsotropicBehaviourIntegrator<
        Hypothesis::PLANESTRAIN>>(ctx, fed, m, std::move(b), params);
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of
     // generatePlaneStrainFaltus2026RegularizedMechanicalBehaviourIntegrators

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<Faltus2026RegularizedIsotropicBehaviourIntegrator<
        Hypothesis::PLANESTRESS>>(ctx, fed, m, std::move(b), params);
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of
     // generatePlaneStressFaltus2026RegularizedMechanicalBehaviourIntegrators

  [[nodiscard]] std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalFaltus2026RegularizedMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    if (b->symmetry != Behaviour::ISOTROPIC) {
      return ctx.registerErrorMessage(
          "only isotropic behaviours are supported");
    }
    auto bi = make_unique<Faltus2026RegularizedIsotropicBehaviourIntegrator<
        Hypothesis::TRIDIMENSIONAL>>(ctx, fed, m, std::move(b), params);
    if (isInvalid(bi)) {
      return {};
    }
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of
     // generateTridimensionalFaltus2026RegularizedMechanicalBehaviourIntegrators

}  // end of namespace mfem_mgis
