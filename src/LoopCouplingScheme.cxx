/*!
 * \file   src/LoopCouplingScheme.cxx
 * \brief  This file implements the inline methods of the `LoopCouplingScheme`
 * class \date   05/12/2022
 */

#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LoopCouplingScheme.hxx"

namespace mfem_mgis {

  std::string LoopCouplingScheme::getDescription() noexcept {
    return {
        "A simple coupling scheme in which each coupling items "
        "is called a fixed number of times"};
  }  // end of  getDescription

  std::map<std::string, std::string>
  LoopCouplingScheme::getParametersDescription() noexcept {
    //     auto d = CouplingSchemeBase::getParametersDescription();
    //     d.insert({"numberOfIterations", "number of iterations of the
    //     scheme"}); return d;
    return {};
  }  // end of getParametersDescription

  LoopCouplingScheme::LoopCouplingScheme(const MeshDiscretization &m)
      : CouplingSchemeBase(m) {}  // end of LoopCouplingScheme

  bool LoopCouplingScheme::setNumberOfIterations(Context &ctx,
                                                 const size_type n) noexcept {
    if (n < 1) {
      return ctx.registerErrorMessage("invalid number of iterations");
    }
    this->number_of_iterations = n;
    return true;
  }  // end of setNumberOfIterations

  //   bool LoopCouplingScheme::addConvergenceCriterion(
  //       Context &ctx, std::string_view, const Parameters &) noexcept {
  //     return ctx.registerErrorMessage(
  //         "the '" + this->getName() +
  //         "' scheme does not support adding a convergence criterion");
  //   }  // end of addConvergenceCriterion

  bool LoopCouplingScheme::addConvergenceCriterion(
      Context &ctx,
      std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>) noexcept {
    return ctx.registerErrorMessage(
        "the '" + this->getName() +
        "' scheme does not support adding a convergence criterion");
  }  // end of addConvergenceCriterion

  std::string LoopCouplingScheme::getName() const noexcept {
    return "LoopCouplingScheme";
  }  // end of getName

  std::optional<std::string> LoopCouplingScheme::describe(
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
    const auto ocis = get_if<bool>(ctx, parameters, "CouplingItems", b);
    if (isInvalid(osd) || isInvalid(onps) || isInvalid(ocis)) {
      return {};
    }
    auto d = std::string{};
    if (*osd) {
      d += "The '" + this->getName() + "' coupling scheme: ";
      d += LoopCouplingScheme::getDescription();
    }
    if (*onps) {
      if (!d.empty()) {
        d += "\n\n";
      }
      d += "The number of iterations per call has been fixed to ";
      d += std::to_string(this->number_of_iterations) + ".";
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

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  LoopCouplingScheme::computeNextState(Context &ctx,
                                       const TimeStep &ts) noexcept {
    if (this->number_of_iterations <= 0) {
      std::ignore = ctx.registerErrorMessage("invalid number of iterations");
      return {ExitStatus::unrecoverableError, {}};
    }
    auto status = ExitStatus{};
    auto iterationsOutputs = std::vector<Parameter>{};
    auto iterationOutputs = std::vector<Parameter>{};
    for (size_type i = 0; i != this->number_of_iterations; ++i) {
      status = ExitStatus{};  // only the last iteration matters
      iterationOutputs.clear();
      ctx.log(verboseLevel2, "* iteration " + std::to_string(i) +
                                 " of the loop coupling scheme");
      for (const auto &m : this->items) {
        ctx.log(verboseLevel2, "* calling computeNextState for '" +
                                   getShortDescription(*m) + "'");
        auto cs = update(ctx, *m);
        const auto o = m->computeNextState(ctx, ts);
        restore(ctx, cs);
        status.update(o.first);
        if (status.shallStop()) {
          ctx.debug("* computeNextState failed for '" +
                    getShortDescription(*m) + "'");
          std::ignore =
              ctx.registerErrorMessage("computing the next state of model '" +
                                       getShortDescription(*m) + "' failed");
          return {status, {}};
        }
        auto itemOutput = Parameters{};
        itemOutput.replaceOrInsert("Name", m->getName());
        itemOutput.replaceOrInsert("Description", getShortDescription(*m));
        if (isValid(o.second)) {
          itemOutput.replaceOrInsert("Output", *(o.second));
        }
        iterationOutputs.push_back(itemOutput);
      }
      iterationsOutputs.push_back(iterationOutputs);
    }
    auto output = ComputeNextStateOutput{};
    output.replaceOrInsert("NumberOfIterations", this->number_of_iterations);
    output.replaceOrInsert("IterationsOutputs", iterationsOutputs);
    output.replaceOrInsert("ItemsOutputs", iterationOutputs);
    //     if ((this->printResourcesUsage) ||  //
    //         (ctx.getVerbosityLevel() >= VerbosityLevel::verboseDebug)) {
    //       print(reportResourcesUsage(profiler, false));
    //     }
    return {status, output};
  }  // end of computeNextState

  LoopCouplingScheme::~LoopCouplingScheme() noexcept = default;

}  // namespace mfem_mgis
