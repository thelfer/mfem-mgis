/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/ResidualOperator.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  NonLinearEvolutionProblem::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : NonLinearEvolutionProblemBase(fed),
        mgis_integrator(new MultiMaterialNonLinearIntegrator(fed, h)),
        hypothesis(h) {
    if (this->fe_discretization->getMesh().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem::revert() {
    NonLinearEvolutionProblemBase::revert();
    this->mgis_integrator->revert();
  }  // end of revert

  void NonLinearEvolutionProblem::update() {
    NonLinearEvolutionProblemBase::update();
    this->mgis_integrator->update();
  }  // end of update

  void NonLinearEvolutionProblem::addBehaviourIntegrator(const std::string& n,
                                                         const size_type m,
                                                         const std::string& l,
                                                         const std::string& b) {
    this->mgis_integrator->addBehaviourIntegrator(n, m, l, b);
  }  // end of addBehaviourIntegrator

  const Material& NonLinearEvolutionProblem::getMaterial(
      const size_type m) const {
    return this->mgis_integrator->getMaterial(m);
  }  // end of NonLinearEvolutionProblem::getMaterial

  Material& NonLinearEvolutionProblem::getMaterial(const size_type m) {
    return this->mgis_integrator->getMaterial(m);
  }  // end of NonLinearEvolutionProblem::getMaterial

  void NonLinearEvolutionProblem::setTimeIncrement(const real dt) {
    this->mgis_integrator->setTimeIncrement(dt);
  }  // end of setTimeIncrement

  NonLinearEvolutionProblem::~NonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis
