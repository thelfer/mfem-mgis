/*!
 * \file   src/ModelBase.cxx
 * \brief  This file implements the `ModelBase` class
 * \date   15/11/2022
 */

#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/ModelBase.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  ModelBase::getParametersDescription() noexcept {
    return getCouplingItemParametersDescription();
  }  // end of getParametersDescription

  ModelBase::ModelBase(PhysicalSystem &ps, const Parameters & /* parameters*/)
      : physicalSystem(ps) {}  // end of ModelBase

  bool ModelBase::completeConstruction(Context &ctx,
                                       const Parameters &parameters) noexcept {
    return handleCouplingItemParameters(ctx, *this, parameters);
  }  // end of completeConstruction

  std::string ModelBase::getIdentifier() const noexcept {
    return this->getName();
  }  // end of getIdentifier

  std::optional<std::string> ModelBase::describe(
      Context &ctx, const bool b, const Parameters &parameters) const noexcept {
    const auto options = std::map<std::string, std::string>{
        {"HelpOptions", "return the options of the describe method"},  //
        {"ShortDescription", "short description"},
        {"DetailedDescription", "detailed description"},
        {"UnknownFields", "list of unknown fields"},
        {"StateVariables", "list of state variables"},
        {"Dependencies", "list of dependencies"}};
    if (!checkParameters(ctx, parameters, options)) {
      return {};
    }
    const auto oho = get_if<bool>(ctx, parameters, "HelpOptions", false);
    if (isInvalid(oho)) {
      return {};
    }
    if (*oho) {
      auto d = std::string{};
      for (auto first = true; const auto &[k, o] : options) {
        if (k == "HelpOptions") {
          continue;
        }
        if (!first) {
          d += '\n';
        }
        d += "- " + k + ": " + o;
        first = false;
      }
      return d;
    }
    const auto osd = get_if<bool>(ctx, parameters, "ShortDescription", true);
    const auto odd = get_if<bool>(ctx, parameters, "DetailedDescription", b);
    const auto oufs = get_if<bool>(ctx, parameters, "UnknownFields", b);
    const auto osvs = get_if<bool>(ctx, parameters, "StateVariables", b);
    const auto odeps = get_if<bool>(ctx, parameters, "Dependencies", b);
    if (areInvalid(osd, odd, oufs, odeps)) {
      return {};
    }
    auto d = std::string{};
    auto add = [&d](const bool ab, const std::string &o) {
      if ((!ab) || (o.empty())) {
        return;
      }
      if (!d.empty()) {
        d += "\n\n";
      }
      d += o;
    };
    add(*osd, "The " + this->getName() + " model");
    add(*odd, this->getDetailedDescription());
    add(*oufs, this->getUnknownFieldsDescription());
    add(*osvs, this->getStateVariablesDescription());
    add(*odeps, this->getDependenciesDescription());
    return d;
  }  // end of describe

  std::string ModelBase::getDetailedDescription() const noexcept {
    return "";
  }  // end of getDetailedDescription

  std::string ModelBase::getUnknownFieldsDescription() const noexcept {
    return "";
  }  // end of getUnknownFieldsDescription

  std::string ModelBase::getStateVariablesDescription() const noexcept {
    auto d = std::string{"This model has no state variable."};
    return d;
  }  // end of getStateVariablesDescription

  std::string ModelBase::getDependenciesDescription() const noexcept {
    return "";
    //     Context ctx;
    //     auto dm = DependenciesManager{};
    //     if (!this->declareDependencies(ctx, dm)) {
    //       return "";
    //     }
    //     const auto a = dm.analyseDependencies();
    //     const auto [ndeps, deps] = getDescription(a);
    //     auto d = std::string{};
    //     if (deps.empty()) {
    //       d += "This model has no mandatory dependency.";
    //     } else if (ndeps == 1) {
    //       d += "This model has the following mandatory dependency:\n";
    //     } else {
    //       d += "This model has the following mandatory dependencies:\n";
    //     }
    //     d += deps;
    //     d += "\n\nNote that this model may declare additional dependencies
    //     during
    //     "; d += "the dependencies resolution process."; return d;
  }  // end of getDependenciesDescription

  VerbosityLevel ModelBase::getVerbosityLevel() const noexcept {
    if (this->verbosityLevel.has_value()) {
      return *(this->verbosityLevel);
    }
    return mgis::getDefaultVerbosityLevel();
  }  // end of getVerbosityLevel

  void ModelBase::setVerbosityLevel(const VerbosityLevel l) noexcept {
    this->verbosityLevel = l;
  }  // end of setVerbosityLevel

  void ModelBase::setLogStream(std::shared_ptr<std::ostream> l) noexcept {
    this->logStream = l;
  }  // end of setLogStream

  std::shared_ptr<std::ostream> ModelBase::getLogStreamPointer() noexcept {
    return this->logStream;
  }  // end of getLogStreamPointer

  PhysicalSystem &ModelBase::getPhysicalSystem() {
    return this->physicalSystem;
  }  // end of getPhysicalSystem

  const PhysicalSystem &ModelBase::getPhysicalSystem() const {
    return this->physicalSystem;
  }  // end of getPhysicalSystem

  std::vector<std::string> ModelBase::getLocations() const noexcept {
    return std::vector<std::string>{};
  }  // end of getLocations

  bool ModelBase::executeInitialPostProcessingTasks(Context &) noexcept {
    return true;
  }  // end of executeInitialPostProcessingTasks

  bool ModelBase::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &) noexcept {
    return true;
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  std::optional<real> ModelBase::getNextTimeIncrement(
      Context &ctx, const real t, const real te) const noexcept {
    auto odt = [this, &ctx, t, te]() -> std::optional<real> {
      auto dt = te - t;
      //     for (const auto &pb : this->bricks) {
      //       ctx.log(verboseLevel3, "** calling `getNextTimeIncrement` on
      //       brick
      //       '", pb->getName(), "'"); const auto odt2 =
      //       pb->getNextTimeIncrement(ctx, t, te); if (isInvalid(odt2)) {
      //         ctx.debug("* computing the next time increment failed for brck
      //         '"
      //         + pb->getName() + "'"); return {};
      //       }
      //       dt = std::min(dt, *odt2);
      //     }
      return dt;
    }();
    //
#pragma message("HERE")
    //    auto b = isTrueOnAllProcesses(isValid(odt));
    auto b = isValid(odt);
    if (!b) {
      return {};
    }
    auto dt = *odt;
#ifdef MFEM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, mpi_type<real>, MPI_MIN,
                  MPI_COMM_WORLD);
#endif /*MFEM_USE_MPI*/
    //
    return dt;
  }  // end of getNextTimeIncrement

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  ModelBase::computeNextState(Context &) noexcept {
    return {ExitStatus::success, ComputeNextStateOutput{}};
  }  // end of computeNextState

  bool ModelBase::executePostProcessingTasks(Context &ctx,
                                             const bool b) noexcept {
    auto r = true;
    for (auto &p : this->postProcessings) {
      if (!p(ctx, b)) {
        r = false;
      }
    }
    return r;
  }  // end of executePostProcessingTasks

  bool ModelBase::update(Context &) noexcept { return true; }  // end update

  bool ModelBase::revert(Context &) noexcept { return true; }  // end revert

  std::vector<std::string> ModelBase::getAvailablePostProcessings()
      const noexcept {
    return {};
  }  // end of getAvailablePostProcessings

  bool ModelBase::addPostProcessing(Context &ctx,
                                    std::string_view n,
                                    const Parameters &) noexcept {
    const auto post_processings = this->getAvailablePostProcessings();
    auto msg = "no post-processing named '" + std::string{n} + "' available";
    msg += ". Avalailable post-processings are:";
    if (!post_processings.empty()) {
      for (const auto &p : post_processings) {
        msg += "\n- " + p;
      }
    }
    return ctx.registerErrorMessage(msg);
  }

  void ModelBase::addPostProcessing(
      std::function<bool(Context &, bool)> p) noexcept {
    this->postProcessings.push_back(p);
  }

  ModelBase::~ModelBase() noexcept = default;

}  // end of namespace mfem_mgis
