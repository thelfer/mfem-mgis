/*!
 * \file   MFEMMGIS/NonLinearModel.cxx
 * \brief
 * \author Thomas Helfer
 * \date   05/03/2026
 */

#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/NonLinearModel.hxx"

namespace mfem_mgis {

  NonLinearModel::NonLinearModel(MeshDiscretization &m,
                                 const Parameters &parameters)
      : ModelBase(m, extract(parameters, ModelBase::getParametersList())),
        problem(std::make_shared<NonLinearEvolutionProblem>(
            m, remove(parameters, ModelBase::getParametersList()))) {
    auto valid_parameters = NonLinearEvolutionProblem::getParametersList();
    for (const auto &[k, d] : ModelBase::getParametersList()) {
      static_cast<void>(d);
      valid_parameters.push_back(k);
    }
    checkParameters(parameters, valid_parameters);
  }

  NonLinearModel::NonLinearModel(std::shared_ptr<NonLinearEvolutionProblem> p)
      : ModelBase(p->getFiniteElementDiscretization()), problem(p) {
    if (p.get() == nullptr) {
      raise("invalid problem");
    }
  }  // end of NonLinearModel

  NonLinearEvolutionProblem &NonLinearModel::getProblem() noexcept {
    return *(this->problem);
  }  // end of getProblem

  const NonLinearEvolutionProblem &NonLinearModel::getProblem() const noexcept {
    return *(this->problem);
  }  // end of getProblem

  // Ajout fonction
  void NonLinearModel::setName(const std::string& n) noexcept {
    this->name = std::move(n);
  }  // end of setName

  // Implémentation fonction
  std::string NonLinearModel::getName() const noexcept {
    return this->name.empty() ? "NonLinearModel" : this->name;
  }  // end of getName

  bool NonLinearModel::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &ctx, const TimeStep &ts) noexcept {
    if (!ModelBase::performInitializationTaksAtTheBeginningOfTheTimeStep(ctx,
                                                                         ts)) {
      return false;
    }
    if (!this->problem->setup(ctx, ts.begin, ts.dt)) {
      return false;
    }
    return true;
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  bool NonLinearModel::executePostProcessingTasks(Context &ctx,
                                                  const TimeStep &ts,
                                                  const bool b) noexcept {
    if (!ModelBase::executePostProcessingTasks(ctx, ts, b)) {
      return false;
    }
    if (b) {
      this->problem->executePostProcessings(ts.begin, ts.end);
    }
    return true;
  }

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  NonLinearModel::computeNextState(Context &ctx, const TimeStep &ts) noexcept {
    const auto r = ModelBase::computeNextState(ctx, ts);
    if (r.first.shallStop()) {
      return {r.first, {}};
    }
    const auto r2 = this->problem->solve(ts.begin, ts.dt);
    if (isInvalid(r2)) {
      return {ExitStatus::recoverableError, {}};
    }
    return {ExitStatus::success, convertToComputeNextStateOutput(r2)};
  }  // end of NonLinearModel

  bool NonLinearModel::update(Context &ctx) noexcept {
    if (!ModelBase::update(ctx)) {
      return false;
    }
    this->problem->update();
    return true;
  }  // end of update

  bool NonLinearModel::revert(Context &ctx) noexcept {
    if (!ModelBase::revert(ctx)) {
      return false;
    }
    this->problem->revert();
    return true;
  }  // end of revert

  NonLinearModel::~NonLinearModel() noexcept = default;

}  // end of namespace mfem_mgis
