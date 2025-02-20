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

  /*!
   * \return a factory for the given hypothesis with some predeclared
   * generators.
   * \tparam H: modelling hypothesis
   *
   * \note this function is meant to be overriden to declare behaviour
   * integrators which are only valid for this modelling hypothesis. Behaviour
   * integrators defined for all modelling hypotheses are declared by the
   * `fillWithDefaultBehaviourIntegrators`.
   */
  template <Hypothesis H>
  static BehaviourIntegratorFactory buildFactory();

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  BehaviourIntegratorFactory buildFactory<Hypothesis::TRIDIMENSIONAL>() {
    BehaviourIntegratorFactory f;
    fillWithDefaultBehaviourIntegrators<Hypothesis::TRIDIMENSIONAL>(f);
    f.addGenerator(
        "Mechanics",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          auto bi = [&]() -> std::unique_ptr<BehaviourIntegrator> {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicTridimensionalStandardFiniteStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          }();
          const auto F = std::array<real, 9u>{1, 1, 1, 0, 0, 0, 0, 0, 0};
          bi->getMaterial().setMacroscopicGradients(F);
          return bi;
        });
    f.addGenerator(
        "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
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
        "TransientHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->btype != Behaviour::GENERALBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          return std::make_unique<TransientHeatTransferBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    f.addGenerator("MicromorphicDamage",  //
                   [](const FiniteElementDiscretization& fed, const size_type m,
                      std::unique_ptr<const Behaviour> b) {
                     return std::make_unique<
                         TridimensionalMicromorphicDamageBehaviourIntegrator>(
                         fed, m, std::move(b));
                   });
    return f;
  }  // end of buildFactory

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  BehaviourIntegratorFactory buildFactory<Hypothesis::PLANESTRAIN>() {
    BehaviourIntegratorFactory f;
    fillWithDefaultBehaviourIntegrators<Hypothesis::PLANESTRAIN>(f);
    f.addGenerator(
        "Mechanics",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicPlaneStrainStandardSmallStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          auto bi = [&]() -> std::unique_ptr<BehaviourIntegrator> {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicPlaneStrainStandardFiniteStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          }();
          const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
          bi->getMaterial().setMacroscopicGradients(F);
          return bi;
        });
    f.addGenerator(
        "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
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
        "MicromorphicDamage",  //
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                BidimensionalMicromorphicDamageBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    return f;
  }  // end of buildFactory

  /*!
   * \brief partial specialisation for the tridimensional case
   */
  template <>
  BehaviourIntegratorFactory buildFactory<Hypothesis::PLANESTRESS>() {
    BehaviourIntegratorFactory f;
    fillWithDefaultBehaviourIntegrators<Hypothesis::PLANESTRESS>(f);
    f.addGenerator(
        "Mechanics",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicPlaneStressStandardSmallStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
            raise("invalid behaviour type");
          }
          auto bi = [&]() -> std::unique_ptr<BehaviourIntegrator> {
            if (b->symmetry == Behaviour::ISOTROPIC) {
              return std::make_unique<
                  IsotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator>(
                  fed, m, std::move(b));
            }
            return std::make_unique<
                OrthotropicPlaneStressStandardFiniteStrainMechanicsBehaviourIntegrator>(
                fed, m, std::move(b));
          }();
          const auto F = std::array<real, 5u>{1, 1, 1, 0, 0};
          bi->getMaterial().setMacroscopicGradients(F);
          return bi;
        });
    f.addGenerator(
        "StationaryNonLinearHeatTransfer",
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
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
        "MicromorphicDamage",  //
        [](const FiniteElementDiscretization& fed, const size_type m,
           std::unique_ptr<const Behaviour> b)
            -> std::unique_ptr<BehaviourIntegrator> {
          if (b->symmetry == Behaviour::ISOTROPIC) {
            return std::make_unique<
                BidimensionalMicromorphicDamageBehaviourIntegrator>(
                fed, m, std::move(b));
          }
          return std::make_unique<
              OrthotropicBidimensionalMicromorphicDamageBehaviourIntegrator>(
              fed, m, std::move(b));
        });
    return f;
  }  // end of buildFactory

  /*!
   * \brief an helper function which add the behaviour integrators for the given
   * hypothesis
   * \tparam H: modelling hypothesis
   * \param[in, out] factories: list of factories per modelling hypotheses
   */
  template <Hypothesis H>
  static void addFactory(
      std::map<Hypothesis, BehaviourIntegratorFactory>& factories) {
    factories.insert({H, buildFactory<H>()});
  }  // end of addFactory

  std::map<Hypothesis, BehaviourIntegratorFactory>
  BehaviourIntegratorFactory::buildFactories() {
    auto factories = std::map<Hypothesis, BehaviourIntegratorFactory>{};
    //     addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN>(factories);
    //     addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS>(factories);
    //     addFactory<Hypothesis::AXISYMMETRICAL>(factories);
    addFactory<Hypothesis::PLANESTRESS>(factories);
    addFactory<Hypothesis::PLANESTRAIN>(factories);
    //     addFactory<Hypothesis::GENERALISEDPLANESTRAIN>(factories);
    addFactory<Hypothesis::TRIDIMENSIONAL>(factories);
    return factories;
  }  // end of BehaviourIntegratorFactory::buildFactories()

  BehaviourIntegratorFactory& BehaviourIntegratorFactory::get(
      const Hypothesis h) {
    static auto factories = BehaviourIntegratorFactory::buildFactories();
    const auto p = factories.find(h);
    if (p == factories.end()) {
      raise("BehaviourIntegratorFactory::get: unsupported hypothesis '" +
            std::string(mgis::behaviour::toString(h)) + "'");
    }
    return p->second;
  }  // end of BehaviourIntegratorFactory::get

  void BehaviourIntegratorFactory::addGenerator(const std::string& n,
                                                const Generator g) {
    if (this->generators.count(n) != 0) {
      raise(
          "BehaviourIntegratorFactory::addGenerator: "
          "a generator named '" +
          n + "' has already been declared");
    }
    this->generators.insert({n, g});
  }  // end of BehaviourIntegratorFactory::addGenerator

  std::unique_ptr<BehaviourIntegrator> BehaviourIntegratorFactory::generate(
      const std::string& n,
      const FiniteElementDiscretization& fed,
      const size_type m,
      std::unique_ptr<const Behaviour> b) const {
    const auto p = this->generators.find(n);
    if (p == this->generators.end()) {
      raise(
          "BehaviourIntegratorFactory::generate: "
          "no generator named '" +
          n + "' declared");
    }
    const auto& g = p->second;
    return g(fed, m, std::move(b));
  }  // end of BehaviourIntegratorFactory::generate

  BehaviourIntegratorFactory::BehaviourIntegratorFactory() = default;
  BehaviourIntegratorFactory::BehaviourIntegratorFactory(
      BehaviourIntegratorFactory&&) = default;
  BehaviourIntegratorFactory::BehaviourIntegratorFactory(
      const BehaviourIntegratorFactory&) = default;
  BehaviourIntegratorFactory& BehaviourIntegratorFactory::operator=(
      BehaviourIntegratorFactory&&) = default;
  BehaviourIntegratorFactory& BehaviourIntegratorFactory::operator=(
      const BehaviourIntegratorFactory&) = default;
  BehaviourIntegratorFactory::~BehaviourIntegratorFactory() noexcept = default;

}  // end of namespace mfem_mgis
