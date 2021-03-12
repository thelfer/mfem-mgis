/*!
 * \file   MultiMaterialEvolutionProblemBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/MultiMaterialEvolutionProblemBase.hxx"

namespace mfem_mgis {

  MultiMaterialEvolutionProblemBase::MultiMaterialEvolutionProblemBase(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : mgis_integrator(new MultiMaterialNonLinearIntegrator(fed, h)),
        hypothesis(h) {}  // end of MultiMaterialEvolutionProblemBase

  void MultiMaterialEvolutionProblemBase::setup(const real t, const real dt) {
    this->mgis_integrator->setup(t, dt);
  }  // end of setup

  void MultiMaterialEvolutionProblemBase::revert() {
    this->mgis_integrator->revert();
  }  // end of revert

  void MultiMaterialEvolutionProblemBase::update() {
    this->mgis_integrator->update();
  }  // end of update

  std::vector<size_type>
  MultiMaterialEvolutionProblemBase::getMaterialIdentifiers() const {
    return this->mgis_integrator->getMaterialIdentifiers();
  }  // end of getMaterialIdentifiers

  void MultiMaterialEvolutionProblemBase::addBehaviourIntegrator(
      const std::string& n,
      const size_type m,
      const std::string& l,
      const std::string& b) {
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& MultiMaterialEvolutionProblemBase::getMaterial(
      const size_type m) const {
    return this->mgis_integrator->getMaterial(m);
  }  // end of MultiMaterialEvolutionProblemBase::getMaterial

  Material& MultiMaterialEvolutionProblemBase::getMaterial(const size_type m) {
    return this->mgis_integrator->getMaterial(m);
  }  // end of MultiMaterialEvolutionProblemBase::getMaterial

  void MultiMaterialEvolutionProblemBase::setTimeIncrement(const real dt) {
    this->mgis_integrator->setTimeIncrement(dt);
  }  // end of setTimeIncrement

  void MultiMaterialEvolutionProblemBase::setMacroscopicGradients(
      mgis::span<const real> g) {
    this->mgis_integrator->setMacroscopicGradients(g);
  }  // end of setMacroscopicGradients

  MultiMaterialEvolutionProblemBase::~MultiMaterialEvolutionProblemBase() =
      default;

}  // end of namespace mfem_mgis
