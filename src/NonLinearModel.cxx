/*!
 * \file   src/NonLinearModel.cxx
 * \brief  This file implements the `NonLinearModel` class
 * \author Thomas Helfer
 * \date   05/03/2026
 */

#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/NonLinearModel.hxx"

namespace mfem_mgis {

  NonLinearModel::NonLinearModel(Context &ctx,
                                 MeshDiscretization &m,
                                 const Parameters &parameters)
      : ModelBase(
            ctx,
            m,
            extract(
                throwing, parameters, ModelBase::getParametersDescription())),
        problem(std::make_shared<NonLinearEvolutionProblem>(
            ctx,
            m,
            remove(parameters, ModelBase::getParametersDescription()))) {
    auto valid_parameters = NonLinearEvolutionProblem::getParametersList();
    for (const auto &[k, d] : ModelBase::getParametersDescription()) {
      static_cast<void>(d);
      valid_parameters.push_back(k);
    }
    checkParameters(throwing, parameters, valid_parameters);
  }

  /*!
   * \return the given problem
   * \param[in] p: non linear evolution problem
   * \throws std::runtime_error if the problem is null
   */
  static NonLinearEvolutionProblem &checkProblem(
      const std::shared_ptr<NonLinearEvolutionProblem> &p) {
    if (p.get() == nullptr) {
      raise("invalid problem");
    }
    return *p;
  }  // end of checkProblem

  NonLinearModel::NonLinearModel(Context &ctx,
                                 std::shared_ptr<NonLinearEvolutionProblem> p)
      : ModelBase(ctx, checkProblem(p).getFiniteElementDiscretization()),
        problem(p) {}  // end of NonLinearModel

  NonLinearEvolutionProblem &NonLinearModel::getProblem() noexcept {
    return *(this->problem);
  }  // end of getProblem

  const NonLinearEvolutionProblem &NonLinearModel::getProblem() const noexcept {
    return *(this->problem);
  }  // end of getProblem

  std::string NonLinearModel::getName() const noexcept {
    return this->name.value_or("NonLinearModel");
  }  // end of getName

  bool NonLinearModel::executeInitialPostProcessingTasks(
      Context &ctx, const real t) noexcept {
    if (!ModelBase::executeInitialPostProcessingTasks(ctx, t)) {
      return false;
    }
    if (!this->problem->executeInitialPostProcessings(ctx, t)) {
      return false;
    }
    return true;
  }  // end of executeInitialPostProcessingTasks

  bool NonLinearModel::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &ctx, const TimeStep &ts) noexcept {
    if (!ModelBase::performInitializationTaksAtTheBeginningOfTheTimeStep(ctx,
                                                                         ts)) {
      return false;
    }
    return true;
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  bool NonLinearModel::executePostProcessingTasks(Context &ctx,
                                                  const TimeStep &ts,
                                                  const bool b) noexcept {
    auto success = ModelBase::executePostProcessingTasks(ctx, ts, b);
    if (b) {
      if (!this->problem->executePostProcessings(ctx, ts.begin, ts.dt)) {
        success = false;
      }
    }
    return success;
  }

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  NonLinearModel::computeNextState(Context &ctx, const TimeStep &ts) noexcept {
    const auto r = ModelBase::computeNextState(ctx, ts);
    if (r.first.shallStop()) {
      return {r.first, {}};
    }
    const auto r2 = this->problem->solve(ctx, ts.begin, ts.dt);
    if (isInvalid(r2)) {
      return {ExitStatus::recoverableError, {}};
    }
    return {ExitStatus::success, convertToComputeNextStateOutput(r2)};
  }  // end of NonLinearModel

  bool NonLinearModel::update(Context &ctx) noexcept {
    if (!ModelBase::update(ctx)) {
      return false;
    }
    return this->problem->update(ctx);
  }  // end of update

  bool NonLinearModel::revert(Context &ctx) noexcept {
    if (!ModelBase::revert(ctx)) {
      return false;
    }
    return this->problem->revert(ctx);
  }  // end of revert

  NonLinearModel::~NonLinearModel() noexcept = default;

}  // end of namespace mfem_mgis
