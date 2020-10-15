/*!
 * \file   src/BehaviourIntegratorFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BehaviourIntegratorFactory.hxx"
#include "MFEMMGIS/SmallStrainMechanicalBehaviourIntegrator.hxx"

namespace mfem_mgis {

  /*!
   * \brief declare behaviour integrators valid for all modelling hypotheses
   *
   * \tparam H: modelling hypothesis
   * \param[in,out] f: factory to be filled
   */
  template <Hypothesis H>
  static void fillWithDefaultBehaviourIntegrators(
      BehaviourIntegratorFactory& f) {
    f.addGenerator(
        "SmallStrainMechanicalBehaviourIntegrator",
        [](const mfem::FiniteElementSpace& fs, const size_type m,
           std::unique_ptr<const Behaviour> b) {
          return std::unique_ptr<SmallStrainMechanicalBehaviourIntegrator<H>>(
              new SmallStrainMechanicalBehaviourIntegrator<H>(fs, m,
                                                              std::move(b)));
        });
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
  static BehaviourIntegratorFactory buildFactory() {
    BehaviourIntegratorFactory f;
    fillWithDefaultBehaviourIntegrators<H>(f);
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
    addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN>(factories);
    addFactory<Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS>(factories);
    addFactory<Hypothesis::AXISYMMETRICAL>(factories);
    addFactory<Hypothesis::PLANESTRESS>(factories);
    addFactory<Hypothesis::PLANESTRAIN>(factories);
    addFactory<Hypothesis::GENERALISEDPLANESTRAIN>(factories);
    addFactory<Hypothesis::TRIDIMENSIONAL>(factories);
    return factories;
  }  // end of BehaviourIntegratorFactory::buildFactories()

  BehaviourIntegratorFactory& BehaviourIntegratorFactory::get(
      const Hypothesis h) {
    static auto factories = BehaviourIntegratorFactory::buildFactories();
    const auto p = factories.find(h);
    if (p == factories.end()) {
      mgis::raise("BehaviourIntegratorFactory::get: unsupported hypothesis '" +
                  std::string(mgis::behaviour::toString(h)) + "'");
    }
    return p->second;
  }  // end of BehaviourIntegratorFactory::get

  void BehaviourIntegratorFactory::addGenerator(const std::string& n,
                                                const Generator g) {
    if (this->generators.count(n) != 0) {
      mgis::raise(
          "BehaviourIntegratorFactory::addGenerator: "
          "a generator named '" +
          n + "' has already been declared");
    }
    this->generators.insert({n, g});
  }  // end of BehaviourIntegratorFactory::addGenerator

  std::unique_ptr<BehaviourIntegrator> BehaviourIntegratorFactory::generate(
      const std::string& n,
      const mfem::FiniteElementSpace& fs,
      const size_type m,
      std::unique_ptr<const Behaviour> b) const {
    const auto p = this->generators.find(n);
    if (p == this->generators.end()) {
      mgis::raise(
          "BehaviourIntegratorFactory::generate: "
          "no generator named '" +
          n + "' declared");
    }
    const auto& g = p->second;
    return g(fs, m, std::move(b));
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