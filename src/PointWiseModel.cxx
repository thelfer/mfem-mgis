/*!
 * \file   PointWiseModel.cxx
 * \brief  This file implements the `PointWiseModel` class
 * \author Thomas Helfer
 * \date   04/05/2026
 */

#include "MGIS/Model/Model.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PointWiseModel.hxx"

namespace mfem_mgis {

  static std::unique_ptr<mgis::model::Model> loadModel(
      attributes::Throwing, const Parameters& parameters) {
    const auto library = get<std::string>(throwing, parameters, "Library");
    const auto model = get<std::string>(throwing, parameters, "Model");
    const auto hypothesis = mgis::behaviour::fromString(
        get<std::string>(throwing, parameters, "Hypothesis"));
    return std::make_unique<mgis::model::Model>(
        mgis::model::load(library, model, hypothesis));
  }  // end of loadModel

  std::map<std::string, std::string>
  PointWiseModel::getParametersDescription() noexcept {
    auto descriptions = ModelBase::getParametersDescription();
    descriptions.insert(
        {{"Library", "path to the library containing the model"},
         {"Model", "model to be loaded"}});
    return descriptions;
  }  // end of getParametersDescription

  PointWiseModel::PointWiseModel(
      std::shared_ptr<const PartialQuadratureSpace> qspace,
      const Parameters& parameters)
      : ModelBase(qspace->getFiniteElementDiscretization()),
        Material(qspace, loadModel(throwing, parameters)) {
    checkParameters(throwing, parameters,
                    PointWiseModel::getParametersDescription());
  }  // end of PointWiseModel

  Material& PointWiseModel::getMaterial() noexcept { return *this; }

  const Material& PointWiseModel::getMaterial() const noexcept { return *this; }

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  PointWiseModel::computeNextState(Context& ctx, const TimeStep& ts) noexcept {
    auto [s, ooutput] = ModelBase::computeNextState(ctx, ts);
    if (!s.shallContinue()) {
      return {s, ooutput};
    }
    const auto opts = mgis::behaviour::BehaviourIntegrationOptions{
        .integration_type =
            mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR,
        .compute_speed_of_sound = false};  // end of BehaviourIntegrationOptions
    const auto r = mgis::behaviour::integrate(*this, opts, ts.dt);
    if (!((r.exit_status == 1) || (r.exit_status == 0))) {
      if (!r.error_message.empty()) {
        std::ignore = ctx.registerErrorMessage(r.error_message);
      }
      return {ExitStatus::recoverableError, {}};
    }
    if (r.exit_status == 0) {
      s.update(ExitStatus::unreliableResults);
    }
    return {s, ooutput};
  }  // end of computeNextState

  bool PointWiseModel::update(Context& ctx) noexcept {
    if (!ModelBase::update(ctx)) {
      return false;
    }
    mgis::behaviour::update(*this);
    return true;
  }  // end of update

  bool PointWiseModel::revert(Context& ctx) noexcept {
    if (!ModelBase::revert(ctx)) {
      return false;
    }
    mgis::behaviour::revert(*this);
    return true;
  }  // end of revert

  PointWiseModel::~PointWiseModel() noexcept = default;

}  // end of namespace mfem_mgis