/*!
 * \file   src/BehaviourIntegratorFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/BidimensionalMicromorphicDamageBehaviourIntegrator.hxx"
#include "MFEMMGIS/OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator.hxx"
#include "MFEMMGIS/TridimensionalMicromorphicDamageBehaviourIntegrator.hxx"
#include "MFEMMGIS/TransientHeatTransferBehaviourIntegrator.hxx"
#include "MFEMMGIS/Faltus2026RegularizedBehaviourIntegrators.hxx"
#include "MFEMMGIS/FBarBehaviourIntegrators.hxx"

namespace mfem_mgis {

  /*!
   * \brief declare behaviour integrators valid for all modelling hypotheses
   *
   * \tparam H: modelling hypothesis
   * \param[in,out] f: factory to be filled
   */
  template <Hypothesis H>
  static void fillWithDefaultBehaviourIntegrators(BehaviourIntegratorFactory&) {
  }  // end of fillWithDefaultBehaviourIntegrators

  [[nodiscard]] static std::unique_ptr<AbstractBehaviourIntegrator>
  generateTridimensionalMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (!checkParameters(
            ctx, params,
            std::map<std::string, std::string>{
                {"Regularization", "Reguralization method (optional)"}})) {
      return {};
    }
    if (contains(params, "Regularization")) {
      const auto oreg_params = get<Parameters>(ctx, params, "Regularization");
      if (isInvalid(oreg_params)) {
        return {};
      }
      const auto ofa = extractFactoryArgument(ctx, *oreg_params);
      if (isInvalid(ofa)) {
        return {};
      }
      if ((ofa->first != "Faltus2026") && (ofa->first != "FBar")) {
        return ctx.registerErrorMessage(
            "invalid regularisation '" + ofa->first +
            "'. The only valid regularisations are 'FBar' and 'Faltus2026'");
      }
      if (ofa->first != "Faltus2026") {
        return generateTridimensionalFaltus2026RegularizedMechanicalBehaviourIntegrators(
            ctx, fed, m, std::move(b), ofa->second);
      }
      return generateTridimensionalFBarBehaviourIntegrators(
          ctx, fed, m, std::move(b), ofa->second);
    }
    if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    auto bi = [&]() -> std::unique_ptr<AbstractBehaviourIntegrator> {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    }();
    const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generateTridimensionalMechanicalBehaviourIntegrators

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  void buildFactory<Hypothesis::TRIDIMENSIONAL>(BehaviourIntegratorFactory& f) {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    fillWithDefaultBehaviourIntegrators<Hypothesis::TRIDIMENSIONAL>(f);
    f.addGenerator(ctx, "Mechanics",
                   generateTridimensionalMechanicalBehaviourIntegrators) |
        or_die;
    f.addGenerator(
        may_abort, "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->btype != Behaviour::GENERALBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                IsotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicTridimensionalStationaryNonLinearHeatTransferBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    f.addGenerator(
        may_abort, "TransientHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->btype != Behaviour::GENERALBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          return std::make_unique<TransientHeatTransferBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    f.addGenerator(may_abort, "MicromorphicDamage",  //
                   [](const FiniteElementDiscretization& fed, const size_type m,
                      std::unique_ptr<const Behaviour> b) {
                     return std::make_unique<
                         TridimensionalMicromorphicDamageBehaviourIntegrator>(
                         fed, m, std::move(b));
                   });
  }  // end of buildFactory

  [[nodiscard]] static std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStrainMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (!checkParameters(
            ctx, params,
            std::map<std::string, std::string>{
                {"Regularization", "Reguralization method (optional)"}})) {
      return {};
    }
    if (contains(params, "Regularization")) {
      const auto oreg_params = get<Parameters>(ctx, params, "Regularization");
      if (isInvalid(oreg_params)) {
        return {};
      }
      const auto ofa = extractFactoryArgument(ctx, *oreg_params);
      if (isInvalid(ofa)) {
        return {};
      }
      if ((ofa->first != "Faltus2026") && (ofa->first != "FBar")) {
        return ctx.registerErrorMessage(
            "invalid regularisation '" + ofa->first +
            "'. The only valid regularisations are 'FBar' and 'Faltus2026'");
      }
      if (ofa->first != "Faltus2026") {
        return generatePlaneStrainFaltus2026RegularizedMechanicalBehaviourIntegrators(
            ctx, fed, m, std::move(b), ofa->second);
      }
      return generatePlaneStrainFBarBehaviourIntegrators(
          ctx, fed, m, std::move(b), ofa->second);
    }
    if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    auto bi = [&]() -> std::unique_ptr<AbstractBehaviourIntegrator> {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    }();
    const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generatePlaneStrainMechanicalBehaviourIntegrators

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  void buildFactory<Hypothesis::PLANESTRAIN>(BehaviourIntegratorFactory& f) {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    fillWithDefaultBehaviourIntegrators<Hypothesis::PLANESTRAIN>(f);
    f.addGenerator(ctx, "Mechanics",
                   generatePlaneStrainMechanicalBehaviourIntegrators) |
        or_die;
    f.addGenerator(
        may_abort, "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->btype != Behaviour::GENERALBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                IsotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicPlaneStrainStationaryNonLinearHeatTransferBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    f.addGenerator(
        may_abort,
        "MicromorphicDamage",  //
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                BidimensionalMicromorphicDamageBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator>(
              fed, m, std::move(b));
        });
  }  // end of buildFactory

  [[nodiscard]] static std::unique_ptr<AbstractBehaviourIntegrator>
  generatePlaneStressMechanicalBehaviourIntegrators(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) noexcept {
    if (!checkParameters(
            ctx, params,
            std::map<std::string, std::string>{
                {"Regularization", "Reguralization method (optional)"}})) {
      return {};
    }
    if (contains(params, "Regularization")) {
      const auto oreg_params = get<Parameters>(ctx, params, "Regularization");
      if (isInvalid(oreg_params)) {
        return {};
      }
      const auto ofa = extractFactoryArgument(ctx, *oreg_params);
      if (isInvalid(ofa)) {
        return {};
      }
      if ((ofa->first != "Faltus2026") && (ofa->first != "FBar")) {
        return ctx.registerErrorMessage(
            "invalid regularisation '" + ofa->first +
            "'. The only valid regularisations are 'FBar' and 'Faltus2026'");
      }
      if (ofa->first != "Faltus2026") {
        return generatePlaneStressFaltus2026RegularizedMechanicalBehaviourIntegrators(
            ctx, fed, m, std::move(b), ofa->second);
      }
      return generatePlaneStressFBarBehaviourIntegrators(
          ctx, fed, m, std::move(b), ofa->second);
    }
    if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      return ctx.registerErrorMessage("invalid behaviour type");
    }
    auto bi = [&]() -> std::unique_ptr<AbstractBehaviourIntegrator> {
      if (b->symmetry == Behaviour::ISOTROPIC) {
        return make_unique<
            IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator>(
            ctx, fed, m, std::move(b));
      }
      return make_unique<
          OrthotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator>(
          ctx, fed, m, std::move(b));
    }();
    const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
    bi->getMaterial().setMacroscopicGradients(F);
    return bi;
  }  // end of generatePlaneStressMechanicalBehaviourIntegrators

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  void buildFactory<Hypothesis::PLANESTRESS>(BehaviourIntegratorFactory& f) {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    fillWithDefaultBehaviourIntegrators<Hypothesis::PLANESTRESS>(f);
    f.addGenerator(ctx, "Mechanics",
                   generatePlaneStressMechanicalBehaviourIntegrators) |
        or_die;
    f.addGenerator(
        may_abort, "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->btype != Behaviour::GENERALBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                IsotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicPlaneStressStationaryNonLinearHeatTransferBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    f.addGenerator(
        may_abort,
        "MicromorphicDamage",  //
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<AbstractBehaviourIntegrator> {
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                BidimensionalMicromorphicDamageBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator>(
              fed, m, std::move(b));
        });
  }  // end of buildFactory

  template <Hypothesis H>
  void BehaviourIntegratorFactory::addFactory(
      std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>&
          factories) {
    auto f = std::unique_ptr<BehaviourIntegratorFactory>(
        new BehaviourIntegratorFactory);
    buildFactory<H>(*f);
    factories.insert({H, std::move(f)});
  }  // end of addFactory

  std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>
  BehaviourIntegratorFactory::buildFactories() {
    auto factories =
        std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>{};
    //     addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN>(factories);
    //     addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS>(factories);
    //     addFactory<Hypothesis::AXISYMMETRICAL>(factories);
    BehaviourIntegratorFactory::addFactory<Hypothesis::PLANESTRESS>(factories);
    BehaviourIntegratorFactory::addFactory<Hypothesis::PLANESTRAIN>(factories);
    //     addFactory<Hypothesis::GENERALISEDPLANESTRAIN>(factories);
    BehaviourIntegratorFactory::addFactory<Hypothesis::TRIDIMENSIONAL>(
        factories);
    return factories;
  }  // end of BehaviourIntegratorFactory::buildFactories()

  OptionalReference<BehaviourIntegratorFactory> BehaviourIntegratorFactory::get(
      Context& ctx, const Hypothesis h) noexcept {
    static auto factories = BehaviourIntegratorFactory::buildFactories();
    const auto p = factories.find(h);
    if (p == factories.end()) {
      return ctx.registerErrorMessage(
          "BehaviourIntegratorFactory::get: unsupported hypothesis '" +
          std::string(mgis::behaviour::toString(h)) + "'");
    }
    return p->second.get();
  }  // end of BehaviourIntegratorFactory::get

  BehaviourIntegratorFactory& BehaviourIntegratorFactory::get(
      const Hypothesis h) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return BehaviourIntegratorFactory::get(ctx, h) | or_raise;
  }  // end of BehaviourIntegratorFactory::get

  bool BehaviourIntegratorFactory::addGenerator(Context& ctx,
                                                const std::string& n,
                                                const Generator g) noexcept {
    if (this->generators.count(n) != 0) {
      return ctx.registerErrorMessage(
          "BehaviourIntegratorFactory::addGenerator: "
          "a generator named '" +
          n + "' has already been declared");
    }
    if (!g) {
      return ctx.registerErrorMessage(
          "BehaviourIntegratorFactory::addGenerator: "
          "invalid generator provided for behaviour integrator '" +
          n + "'");
    }
    this->generators.insert({n, g});
    return true;
  }  // end of addGenerator

  void BehaviourIntegratorFactory::addGenerator(
      attributes::MayAbort,
      std::string_view n,
      const DeprecatedGeneratorType g) {
    auto ng = [g](Context& nctx, const FiniteElementDiscretization& nfed,
                  const size_type nm, std::unique_ptr<const Behaviour> nb,
                  const Parameters& nparams)
        -> std::unique_ptr<AbstractBehaviourIntegrator> {
      if (!nparams.empty()) {
        return nctx.registerErrorMessage("no parameter expected");
      }
      try {
        return g(nfed, nm, std::move(nb));
      } catch (...) {
        std::ignore = registerExceptionInErrorBacktrace(nctx);
      }
      return {};
    };
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    this->addGenerator(ctx, std::string{n}, ng) | or_die;
  }  // end of addGenerator

  std::unique_ptr<AbstractBehaviourIntegrator>
  BehaviourIntegratorFactory::generate(
      Context& ctx,
      std::string_view n,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b,
      const Parameters& params) const noexcept {
    const auto p = this->generators.find(n);
    if (p == this->generators.end()) {
      return ctx.registerErrorMessage(
          "BehaviourIntegratorFactory::generate: "
          "no generator named '" +
          std::string{n} + "' declared");
    }
    const auto& g = p->second;
    try {
      auto bi = g(ctx, fed, m, std::move(b), params);
      if (isInvalid(bi)) {
        return ctx.registerErrorMessage(
            "BehaviourIntegratorFactory::generate: "
            "generation of behaviour integrator '" +
            std::string{n} + "' failed");
      }
      return bi;
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of BehaviourIntegratorFactory::generate

  std::unique_ptr<AbstractBehaviourIntegrator>
  BehaviourIntegratorFactory::generate(
      std::string_view n,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b) const {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return this->generate(ctx, n, fed, m, std::move(b), {}) | or_raise;
  }  // end of BehaviourIntegratorFactory::generate

  BehaviourIntegratorFactory::BehaviourIntegratorFactory() = default;
  BehaviourIntegratorFactory::BehaviourIntegratorFactory(
      BehaviourIntegratorFactory&&) = default;
  BehaviourIntegratorFactory::~BehaviourIntegratorFactory() noexcept = default;

}  // end of namespace mfem_mgis
