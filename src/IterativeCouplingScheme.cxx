/*!
 * \file   src/IterativeCouplingScheme.cxx
 * \brief  This file implements the inline methods of the
 * `IterativeCouplingScheme` class \date   05/12/2022
 */

#include "MFEMMGIS/AbstractCouplingSchemeConvergenceCriterion.hxx"
#include "MFEMMGIS/IterativeCouplingScheme.hxx"

namespace mfem_mgis {

  std::string IterativeCouplingScheme::getDescription() noexcept {
    return {
        "A coupling scheme which calls each coupling items "
        "until some convergence criteria are satisfied"};
  }  // end of  getDescription

  std::map<std::string, std::string>
  IterativeCouplingScheme::getParametersDescription() noexcept {
    auto d = CouplingSchemeBase::getParametersDescription();
    //     d.insert(
    //         {"MaximumNumberOfIterations", "number of iterations of the
    //         scheme"});
    //     d.insert({"ConvergenceCriteria",
    //               "list of convergence criteria (array of Parameters)"});
    //     return d;
    return {};
  }  // end of getParametersDescription

  IterativeCouplingScheme::IterativeCouplingScheme() = default;

  bool IterativeCouplingScheme::setMaximumNumberOfIterations(
      Context &ctx, const size_type n) noexcept {
    if (n < 1) {
      return ctx.registerErrorMessage("invalid number of iterations");
    }
    this->maximum_number_of_iterations = n;
    return true;
  }  // end of setMaximumNumberOfIterations

  //   bool IterativeCouplingScheme::addConvergenceCriterion(
  //       Context &ctx, std::string_view n, const Parameters &parameters)
  //       noexcept {
  //     auto &f = CouplingSchemeConvergenceCriterionFactory::get();
  //     auto c = f.create(ctx, std::string{n}, this->physicalSystem_,
  //     parameters); return this->addConvergenceCriterion(ctx, c);
  //   }

  bool IterativeCouplingScheme::addConvergenceCriterion(
      Context &ctx,
      std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion> c) noexcept {
    if (c.get() == nullptr) {
      return ctx.registerErrorMessage("invalid criteria");
    }
    this->convergence_criteria.push_back(c);
    return true;
  }  // end of addConvergenceCriterion

  std::string IterativeCouplingScheme::getName() const noexcept {
    return "IterativeCouplingScheme";
  }  // end of getName

  std::optional<std::string> IterativeCouplingScheme::describe(
      Context &ctx, const bool b, const Parameters &parameters) const noexcept {
    if (!checkParameters(
            ctx,
            parameters,  //
            std::map<std::string, std::string>{
                {"ShortDescription", "short description"},
                {"NumericalParameters", "list of numerical parameters"},
                {"CouplingItems", "list of coupling items"}})) {
      return {};
    }
    const auto osd = get_if<bool>(ctx, parameters, "ShortDescription", true);
    const auto onps = get_if<bool>(ctx, parameters, "NumericalParameters", b);
    const auto ocis = get_if<bool>(ctx, parameters, "NumericalParameters", b);
    if (isInvalid(osd) || isInvalid(onps) || isInvalid(ocis)) {
      return {};
    }
    auto d = std::string{};
    if (*osd) {
      d += "The '" + this->getName() + "' coupling scheme: ";
      d += IterativeCouplingScheme::getDescription();
    }
    if (*onps) {
      if (!d.empty()) {
        d += "\n\n";
      }
      d += "The maximum number of iterations per call has been fixed to ";
      d += std::to_string(this->maximum_number_of_iterations) + ".";
    }
    if (*ocis) {
      const auto cid = this->getCouplingItemsDescription();
      if (!cid.empty()) {
        if (!d.empty()) {
          d += "\n\n";
        }
        d += cid;
      }
    }
    return d;
  }  // end of describe

  bool
  IterativeCouplingScheme::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &ctx) noexcept {
    const auto r = CouplingSchemeBase::
        performInitializationTaksAtTheBeginningOfTheTimeStep(ctx);
    if (!r) {
      return false;
    }
    for (const auto &c : this->convergence_criteria) {
      const auto r2 =
          c->performInitializationTaksAtTheBeginningOfTheTimeStep(ctx);
      if (!r2) {
        ctx.debug(
            "* performInitializationTaksAtTheBeginningOfTheTimeStep failed for "
            "convergence criterion");
        return false;
      }
    }
    return true;
  }

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  IterativeCouplingScheme::computeNextState(Context &ctx) noexcept {
    //     [[maybe_unused]] auto profiler =
    //         ctx.startResourcesProfiling<1>("LoopCouplingScheme::computeNextState");
    if (this->convergence_criteria.empty()) {
      std::ignore =
          ctx.registerErrorMessage("no convergence criterion defined");
      return {ExitStatus::unrecoverableError, {}};
    }
    if (this->maximum_number_of_iterations <= 0) {
      std::ignore =
          ctx.registerErrorMessage("invalid maximum number of iterations");
      return {ExitStatus::unrecoverableError, {}};
    }
    auto status = ExitStatus{};
    auto iterationsOutputs = std::vector<Parameter>{};
    for (size_type i = 0; i != this->maximum_number_of_iterations; ++i) {
      status = ExitStatus{};  // only the last iteration matters
      auto iterationOutputs = std::vector<Parameter>{};
      ctx.log(verboseLevel2, "* iteration " + std::to_string(i) +
                                 " of the iterative coupling scheme");
      for (const auto &m : this->items) {
        //         [[maybe_unused]] auto iprofiler =
        //             ctx.startResourcesProfiling<1>(m->getName() +
        //             "::computeNextState");
        ctx.log(verboseLevel2, "* calling computeNextState for '" +
                                   getShortDescription(*m) + "'");
        auto cs = CouplingSchemeBase::update(ctx, *m);
        const auto o = m->computeNextState(ctx);
        restore(ctx, cs);
        status.update(o.first);
        if (status.shallStop()) {
          ctx.debug("* computeNextState failed for '" +
                    getShortDescription(*m) + "'");
          std::ignore =
              ctx.registerErrorMessage("computing the next state of model '" +
                                       getShortDescription(*m) + "' failed");
          return {o.first, {}};
        }
        auto itemOutput = Parameters{};
        itemOutput.replaceOrInsert("Name", m->getName());
        itemOutput.replaceOrInsert("Description", getShortDescription(*m));
        if (isValid(o.second)) {
          itemOutput.replaceOrInsert("Output", *(o.second));
        }
        iterationOutputs.push_back(itemOutput);
      }  // end of loop on coupling items
      iterationsOutputs.push_back(iterationOutputs);
      // preparing the output of the coupling scheme
      auto output = ComputeNextStateOutput{};
      output.replaceOrInsert("NumberOfIterations", i);
      output.replaceOrInsert("IterationOutputs", iterationsOutputs);
      output.replaceOrInsert("ItemsOutputs", iterationOutputs);
      //
      auto converged = true;
      for (const auto &c : this->convergence_criteria) {
        const auto ob = c->check(ctx, output);
        if (isInvalid(ob)) {
          return {ExitStatus::unrecoverableError, {}};
        }
        if (!(*ob)) {
          converged = false;
        }
      }
      if (converged) {
        ctx.log(verboseLevel2, "* iterative coupling scheme '" +
                                   this->getName() + "' converged after " +
                                   std::to_string(i + 1) + " iterations");
        //    if ((this->printResourcesUsage_) ||  //
        //             (ctx.getVerbosityLevel() >=
        //             VerbosityLevel::verboseDebug)) {
        //           print(reportResourcesUsage(profiler, true));
        //         }
        return {status, output};
      }
    }
    ctx.log(verboseLevel2,
            "* iterative coupling scheme '" + this->getName() +
                "'did no converge after " +
                std::to_string(this->maximum_number_of_iterations + 1) +
                " iterations");
    //     if ((this->printResourcesUsage_) ||  //
    //         (ctx.getVerbosityLevel() >= VerbosityLevel::verboseDebug)) {
    //       print(reportResourcesUsage(profiler, false));
    //     }
    return {ExitStatus::recoverableError, {}};
  }  // end of computeNextState

  bool IterativeCouplingScheme::update(Context &ctx) noexcept {
    if (!CouplingSchemeBase::update(ctx)) {
      return false;
    }
    for (const auto &c : this->convergence_criteria) {
      const auto s = c->update(ctx);
      if (!s) {
        ctx.debug("* update failed for convergence criterion");
        return false;
      }
    }
    return true;
  }  // end of update

  bool IterativeCouplingScheme::revert(Context &ctx) noexcept {
    if (!CouplingSchemeBase::revert(ctx)) {
      return false;
    }
    for (const auto &c : this->convergence_criteria) {
      const auto s = c->revert(ctx);
      if (!s) {
        ctx.debug("* revert failed for convergence criterion");
        return false;
      }
    }
    return true;
  }  // end of revert

  IterativeCouplingScheme::~IterativeCouplingScheme() noexcept = default;

}  // namespace mfem_mgis
