/*!
 * \file   src/Simulation.cxx
 * \brief  This file implements the methods of the Simulation class.
 * \date   15/05/2023
 */

#include <cmath>
#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/DefaultTimeStepValidator.hxx"
#include "MFEMMGIS/DefaultTimeIncrementComputer.hxx"
#include "MFEMMGIS/DefaultConvergenceFailureHandler.hxx"
#include "MFEMMGIS/Simulation.hxx"

namespace mfem_mgis {

  void ExitStatus::synchronize() noexcept {
#ifdef MFEM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &(this->status.value), 1, mpi_type<size_type>,
                  MPI_SUM, MPI_COMM_WORLD);
#endif /* MFEM_USE_MPI */
  }    // end of synchronize

  ExitStatus synchronize(const ExitStatus s) noexcept {
    auto ns = s;
    ns.synchronize();
    return ns;
  }  // end of synchronize

  SimulationOutput::SimulationOutput() noexcept = default;
  SimulationOutput::SimulationOutput(SimulationOutput &&) noexcept = default;
  SimulationOutput::SimulationOutput(const SimulationOutput &) noexcept =
      default;
  SimulationOutput &SimulationOutput::operator=(SimulationOutput &&) noexcept =
      default;
  SimulationOutput &SimulationOutput::operator=(
      const SimulationOutput &) noexcept = default;
  SimulationOutput::~SimulationOutput() noexcept = default;

  Simulation::TimesDescription::TimesDescription(const real tb, const real te)
      : TimesDescription(tb, te, 1) {}  // end of TimesDescription

  Simulation::TimesDescription::TimesDescription(const real tb,
                                                 const real te,
                                                 const size_type n) {
    if (tb > te) {
      raise("invalid time sequence between times '" + std::to_string(tb) +
            "' and '" + std::to_string(te) + "'");
    }
    if (n < 1) {
      raise("invalid number of temporal sequences (" + std::to_string(n) + ")");
    }
    const auto dt = (te - tb) / n;
    if (std::fpclassify((tb + dt) - tb) == FP_ZERO) {
      raise("time increment too small (" + std::to_string(dt) + ")");
    }
    this->push_back(tb);
    auto t = tb + dt;
    for (size_type i = 0; i != n - 1; ++i) {
      this->push_back(t);
      t += dt;
    }
    this->push_back(te);
  }  // end of TimesDescription

  Simulation::TimesDescription::TimesDescription(std::span<const real> times) {
    if (times.size() < 2) {
      raise("invalid number of times");
    }
    auto p_bts = times.begin();
    auto p_ets = std::next(p_bts);
    auto pe = times.end();
    while (p_ets != pe) {
      if (*p_ets - *p_bts < 0) {
        raise("negative time increment between times '" +
              std::to_string(*p_bts) + "' and '" + std::to_string(*p_ets) +
              "'");
      }
      if (std::fpclassify(*p_ets - *p_bts) == FP_ZERO) {
        raise("invalid time increment at time '" + std::to_string(*p_bts) +
              "'");
      }
      ++p_bts;
      ++p_ets;
    }
    std::vector<real>::insert(this->end(), times.begin(), times.end());
  }

  std::map<std::string, std::string>
  Simulation::getParametersDescription() noexcept {
    auto d = std::map<std::string, std::string>{};
    d.insert({"Times", "list of times defining the temporal sequences"});
    d.insert({"Monitors", "list of monitors (vector of parameters)"});
    d.insert({"KeepOutputs", "keep the simulation outputs for all time steps"});
    d.insert({"allowSubStepping",
              "boolean stating if sub stepping in case of convergence failure "
              "is allowed"});
    d.insert({"IndependentTemporalSequences",
              "boolean stating if the temporal sequences are independent. "
              "Currently, this boolean only choose if the last time "
              "increment of the previous temporal sequence is taken into "
              "account to bound the first time increment of the "
              "current temporal sequence (bounding occurs if this boolean is "
              "false)"});
    d.insert({"MinimalTimeIncrement", "minimal time increment"});
    d.insert({"MaximalTimeIncrement", "maximal time increment"});
    d.insert({"MaximumNumberOfFailuresPerTemporalSequence",
              "maximum number of failures or time step rejections allowed "
              "within each temporal sequence"});
    d.insert(
        {"TimeStepValidator", "strategy used to determine the next time step"});
    d.insert({"timeIncrementComputer",
              "strategy used to determine the next time step"});
    d.insert(
        {"LimitTimeIncrementIncrease",
         "boolean stating if the current estimate of the next time step can "
         "be greater than the previous time increment "
         "multiplied by the 'maximalTimeIncrementRelativeIncrease' parameter"});
    d.insert({"MaximalTimeIncrementRelativeIncrease",
              "coefficient used to determined the maximum ratio between the "
              "next time step "
              "and the previous one"});
    d.insert({"LimitTimeIncrementDecrease",
              "boolean stating if the current estimate of the next time step "
              "can be lower "
              "than the previous time increment multiplied by the "
              "'maximalTimeIncrementRelativeDecrease' parameter"});
    d.insert({"MaximalTimeIncrementRelativeDecrease",
              "coefficient used to determined the minimum ratio between the "
              "next time step "
              "and the previous one"});
    d.insert({"TimeIncrementBalancer",
              "strategy used to balance the time increments to avoid small "
              "time increments "
              "at the end of the temporal sequence"});
    d.insert({"BalanceTimeIncrement",
              "boolean stating if the time steps must be balanced to avoid "
              "small time increments "
              "at the end of the temporal sequence"});
    d.insert({"ConvergenceFailureHandler",
              "strategy used to determine how a divergence of "
              "the resolution shall be handled"});
    d.insert({"MaximumNumberOfTimeSteps",
              "the maximum number of time steps allowed per call to the run "
              "method"});
    d.insert({"NumberOfTimeStepsBetweenPostProcessings",
              "the number of time steps between two post-processings marked as "
              "'explicitly requested by the user'. Every time a "
              "post-processing is called, it recieves a boolean stating if the "
              "post-processing time as been 'explicitly "
              "requested by the user'. Lightweight post-processings (`Curves` "
              "for instance) may ignore this boolean. This "
              "boolean is mostly important for heavy post-processsings in "
              "terms of memory, disk usage or computations, such "
              "as the `VTKExport` post-processing. Note that the each end of a "
              "temporal sequence is always flaged as 'explicitly "
              "requested by the user', independently of the "
              "`numberOfTimeStepsBetweenPostProcessings` parameter. Note also "
              "that the "
              "number of "
              "computed time steps is reset at each call of the `run` method. "
              "The parameter `numberOfTimeStepsBetweenPostProcessings` "
              "thus has no effect if it is greated than "
              "`maximumNumberOfTimeSteps` parameter."});
    d.insert({"TimeBetweenPostProcessings",
              "the time between two post-processings marked as 'explicitly "
              "requested by the user'. Every time a "
              "post-processing is called, it recieves a boolean stating if the "
              "post-processing time as been 'explicitly "
              "requested by the user'. Lightweight post-processings (`Curves` "
              "for instance) may ignore this boolean. This "
              "boolean is mostly important for heavy post-processsings in "
              "terms of memory, disk usage or computations, such "
              "as the `VTKExport` post-processing. Note that the each end of a "
              "temporal sequence is always flaged as 'explicitly "
              "requested by the user', independently of the "
              "`timeBetweenPostProcessings` parameter. Note also that the time "
              "since the beginning of the simulation is reset at each call of "
              "the `run` method."});
    return d;
  }  // end of getParametersDescription

  Simulation::TimesDescription::TimesDescription(TimesDescription &&) noexcept =
      default;
  Simulation::TimesDescription::TimesDescription(
      const TimesDescription &) noexcept = default;
  Simulation::TimesDescription &Simulation::TimesDescription::operator=(
      TimesDescription &&) noexcept = default;
  Simulation::TimesDescription &Simulation::TimesDescription::operator=(
      const TimesDescription &) noexcept = default;

  Simulation::TimesDescription::~TimesDescription() noexcept = default;

  [[nodiscard]] static std::vector<real> getTimes(
      attributes::Throwing, const Parameters &parameters) {
    auto r = std::vector<real>{};
    const auto times =
        get<std::vector<Parameter>>(throwing, parameters, "times");
    r.reserve(times.size());
    for (const auto &v : times) {
      r.push_back(get<real>(throwing, v));
    }
    return r;
  }

  Simulation::Simulation(AbstractNonLinearEvolutionProblem &p,
                         const Parameters &parameters)
      : nonlinearEvolutionProblem(&p),
        timesDescription(::mfem_mgis::getTimes(throwing, parameters)),
        keepOutputs(get_if(throwing, parameters, "keepOutputs", false)) {
    this->treatParameters(throwing, parameters);
    this->completeInitialization();
  }  // end of Simulation

  Simulation::Simulation(AbstractNonLinearEvolutionProblem &p,
                         const TimesDescription &times)
      : nonlinearEvolutionProblem(&p),
        timesDescription(times),
        keepOutputs(false) {
    this->completeInitialization();
  }  // end of Simulation

  Simulation::Simulation(AbstractNonLinearEvolutionProblem &p,
                         const std::initializer_list<real> &times)
      : nonlinearEvolutionProblem(&p),
        timesDescription(times),
        keepOutputs(false) {
    this->completeInitialization();
  }  // end of Simulation

  // Simulation::Simulation(PhysicalSystem &ps, const Parameters &parameters)
  //     : physicalSystem(ps)
  //     , timesDescription(get<std::vector<real>>(parameters, "times"))
  //     , keepOutputs(get_if(parameters, "keepOutputs", false))
  // {
  //   this->treatParameters(attributes::throwing, params);
  //   this->completeInitialization();
  // }    // end of Simulation
  //
  // Simulation::Simulation(PhysicalSystem &ps, const TimesDescription &times)
  //     : physicalSystem(ps), timesDescription(times), keepOutputs(false)
  // {
  //   this->completeInitialization();
  // }    // end of Simulation
  //
  // Simulation::Simulation(PhysicalSystem &ps, const
  // std::initializer_list<real> &times)
  //     : physicalSystem(ps), timesDescription(times), keepOutputs(false)
  // {
  //   this->completeInitialization();
  // }    // end of Simulation
  //

  void Simulation::treatParameters(attributes::Throwing,
                                   const Parameters &parameters) {
    checkParameters(throwing, parameters,
                    Simulation::getParametersDescription());
    this->allowSubStepping =
        get_if<bool>(throwing, parameters, "allowSubStepping", true);
    this->independentTemporalSequences = get_if<bool>(
        throwing, parameters, "independentTemporalSequences", true);
    if (contains(parameters, "minimalTimeIncrement")) {
      this->minimalTimeIncrement =
          get<real>(throwing, parameters, "minimalTimeIncrement");
      if ((*(this->minimalTimeIncrement) < 0) ||
          (std::fpclassify(*(this->minimalTimeIncrement)) == FP_ZERO)) {
        raise("null or negative minimal time increment (" +
              std::to_string(*(this->minimalTimeIncrement)) + ")");
      }
    }
    if (contains(parameters, "maximalTimeIncrement")) {
      this->maximalTimeIncrement =
          get<real>(throwing, parameters, "maximalTimeIncrement");
      if ((*(this->maximalTimeIncrement) < 0) ||
          (std::fpclassify(*(this->maximalTimeIncrement)) == FP_ZERO)) {
        raise("negative maximal time increment (" +
              std::to_string(*(this->maximalTimeIncrement)) + ")");
      }
    }
    if ((this->minimalTimeIncrement.has_value()) &&
        (this->maximalTimeIncrement.has_value())) {
      if (*(this->minimalTimeIncrement) > *(this->maximalTimeIncrement)) {
        raise("minimal time increment (" +
              std::to_string(*(this->minimalTimeIncrement)) +
              ") is greater than maximal time increment (" +
              std::to_string(*(this->maximalTimeIncrement)) + ")");
      }
    }
    if (this->minimalTimeIncrement.has_value()) {
      const auto pe = this->timesDescription.end();
      auto p_bts = this->timesDescription.begin();
      auto p_ets = std::next(p_bts);
      while (p_ets != pe) {
        if (*p_ets - *p_bts < *(this->minimalTimeIncrement)) {
          raise("time sequence (" + std::to_string(*p_bts) + ", " +
                std::to_string(*p_ets) +
                ") is not compatible with the minimal time increment (" +
                std::to_string(*(this->minimalTimeIncrement)) + ")");
        }
        ++p_bts;
        ++p_ets;
      }
    }
    if (contains(parameters, "maximumNumberOfFailuresPerTemporalSequence")) {
      this->maximumNumberOfFailuresPerTemporalSequence = get<size_type>(
          throwing, parameters, "maximumNumberOfFailuresPerTemporalSequence");
      if (this->maximumNumberOfFailuresPerTemporalSequence < 1) {
        raise("invalid number of failures per temporal sequence (" +
              std::to_string(this->maximumNumberOfFailuresPerTemporalSequence) +
              ")");
      }
    }
    this->limitTimeIncrementIncrease =
        get_if<bool>(throwing, parameters, "limitTimeIncrementIncrease", true);
    if (contains(parameters, "maximalTimeIncrementRelativeIncrease")) {
      if (!this->limitTimeIncrementIncrease) {
        raise(
            "defining the 'maximalTimeIncrementRelativeIncrease' parameter is "
            "meaningless if the 'limitTimeIncrementIncrease' "
            "parameter is set to false");
      }
      this->maximalTimeIncrementRelativeIncrease = get<real>(
          throwing, parameters, "maximalTimeIncrementRelativeIncrease");
      if (this->maximalTimeIncrementRelativeIncrease < 1) {
        raise(
            "the 'maximalTimeIncrementRelativeIncrease' parameter must be "
            "greater than 1");
      }
    }
    this->limitTimeIncrementDecrease =
        get_if<bool>(throwing, parameters, "limitTimeIncrementDecrease", true);
    if (contains(parameters, "maximalTimeIncrementRelativeDecrease")) {
      if (!this->limitTimeIncrementDecrease) {
        raise(
            "defining the 'maximalTimeIncrementRelativeDecrease' parameter is "
            "meaningless if the 'limitTimeIncrementDecrease' "
            "parameter is set to false");
      }
      this->maximalTimeIncrementRelativeDecrease = get<real>(
          throwing, parameters, "maximalTimeIncrementRelativeDecrease");
      if (this->maximalTimeIncrementRelativeDecrease > 1) {
        raise(
            "the 'maximalTimeIncrementRelativeDecrease' parameter must be "
            "lower than 1");
      }
    }
    if (contains(parameters, "balanceTimeIncrement")) {
      this->balanceTimeIncrements =
          get<bool>(throwing, parameters, "balanceTimeIncrement");
    }
    if (contains(parameters, "timeIncrementBalancer")) {
      if (!this->balanceTimeIncrements) {
        raise(
            "specifying a time increment balancer is meaningless if time "
            "increment balancing is disabled");
      }
      const auto [n, tbparams] = extractFactoryArgument(
          throwing,
          get<Parameters>(throwing, parameters, "timeIncrementBalancer"));
      if (n != "Default") {
        raise("unsupported time increment balancer '" + n + "'");
      }
      auto validParameters = std::map<std::string, std::string>{
          {"minimalRelativeRemainder",
           "value used to check if the rest of the time sequence is exactly "
           "the product of the proposed time step by an integer"},
          {"maximalRelativeRemainder",
           "value used to decide if it is worth to balance the time steps"}};
      checkParameters(throwing, tbparams, validParameters);
      if (contains(tbparams, "minimalRelativeRemainder")) {
        this->minimalRelativeRemainder =
            get<real>(throwing, tbparams, "minimalRelativeRemainder");
      }
      if (contains(tbparams, "maximalRelativeRemainder")) {
        this->maximalRelativeRemainder =
            get<real>(throwing, tbparams, "maxmalRelativeRemainder");
      }
      auto checkBounds = [](const real r, const auto *const rn) {
        if (r >= 0.99999) {
          raise("the '" + std::string{rn} +
                "'relative remainder must be lower than 0.99999");
        }
        if (r <= 1.e-6) {
          raise("the '" + std::string{rn} +
                "'relative remainder must be greater than 10e-6");
        }
      };
      checkBounds(this->minimalRelativeRemainder, "minimal");
      checkBounds(this->maximalRelativeRemainder, "maximal");
      if (this->maximalRelativeRemainder < this->minimalRelativeRemainder) {
        raise(
            "maximal relative remainder must be greater than the minimal one");
      }
    }
    if (contains(parameters, "timeStepValidator")) {
      const auto [n, tsv] = extractFactoryArgument(
          throwing, get<Parameters>(throwing, parameters, "timeStepValidator"));
      if (n != "Default") {
        raise("unsupported time validator '" + n + "'");
      }
      if (!tsv.empty()) {
        raise("no parameters expected for time validator '" + n + "'");
      }
    }
    if (contains(parameters, "timeIncrementComputer")) {
      const auto [n, tic] = extractFactoryArgument(
          throwing,
          get<Parameters>(throwing, parameters, "timeIncrementComputer"));
      if (n != "Default") {
        raise("unsupported time increment computer '" + n + "'");
      }
      if (!tic.empty()) {
        raise("no parameters expected for time increment computer '" + n + "'");
      }
    }
    if (contains(parameters, "convergenceFailureHandler")) {
      const auto [n, tic] = extractFactoryArgument(
          throwing,
          get<Parameters>(throwing, parameters, "convergenceFailureHandler"));
      if (n != "Default") {
        raise("unsupported convergence failure handler '" + n + "'");
      }
      if (!tic.empty()) {
        raise("no parameters expected for convergence failure handler '" + n +
              "'");
      }
    }
    if (contains(parameters, "maximumNumberOfTimeSteps")) {
      this->maximumNumberOfTimeSteps =
          get<size_type>(throwing, parameters, "maximumNumberOfTimeSteps");
      if (*(this->maximumNumberOfTimeSteps) < 1) {
        raise("invalid maximum number of time steps (" +
              std::to_string(*(this->maximumNumberOfTimeSteps)) + ")");
      }
    }
    if (contains(parameters, "numberOfTimeStepsBetweenPostProcessings")) {
      this->numberOfTimeStepsBetweenPostProcessings = get<size_type>(
          throwing, parameters, "numberOfTimeStepsBetweenPostProcessings");
      if (*(this->numberOfTimeStepsBetweenPostProcessings) < 1) {
        raise("invalid number of time steps between post-processings (" +
              std::to_string(*(this->numberOfTimeStepsBetweenPostProcessings)) +
              ")");
      }
    }
    if (contains(parameters, "timeBetweenPostProcessings")) {
      this->timeBetweenPostProcessings =
          get<real>(throwing, parameters, "timeBetweenPostProcessings");
      if (*(this->timeBetweenPostProcessings) <
          std::numeric_limits<real>::min()) {
        raise("invalid time between post-processings " +
              std::to_string(*(this->timeBetweenPostProcessings)) + ")");
      }
    }
    //     if (contains(parameters, "monitors")) {
    //       auto ctx = Context{};
    //       for (const auto &m :
    //            get<std::vector<Parameter>>(throwing, parameters, "monitors"))
    //            {
    //         const auto [n, mparams] =
    //             extractFactoryArgument(throwing, get<Parameters>(throwing,
    //             m));
    //         if (!this->addMonitor(ctx, n, mparams)) {
    //           raise(ctx.getErrorMessage());
    //         }
    //       }
    //     }
  }  // end of treatParameters

  void Simulation::completeInitialization() {
    if (this->timeStepValidator == nullptr) {
      this->timeStepValidator = std::make_shared<DefaultTimeStepValidator>();
    }
    if (this->timeIncrementComputers.empty()) {
      auto computer = std::make_shared<DefaultTimeIncrementComputer>();
      this->timeIncrementComputers.push_back(computer);
    }
    if (this->convergenceFailureHandler == nullptr) {
      this->convergenceFailureHandler =
          std::make_shared<DefaultConvergenceFailureHandler>();
    }
  }  // end of completeInitialization

  void Simulation::addInitializationTask(const InitializationTask &t) noexcept {
    this->addInitializationTask(
        std::to_string(this->initializationTasks.size()), t);
  }  // end of addPostProcessingTask

  void Simulation::addInitializationTask(const std::string &n,
                                         const InitializationTask &t) noexcept {
    this->initializationTasks.push_back({n, t});
  }  // end of addPostProcessingTask

  void Simulation::addPostProcessingTask(const PostProcessingTask &t) noexcept {
    this->addPostProcessingTask(
        std::to_string(this->postProcessingTasks.size()), t);
  }  // end of addPostProcessingTask

  void Simulation::addPostProcessingTask(const std::string &n,
                                         const PostProcessingTask &t) noexcept {
    this->postProcessingTasks.push_back({n, t});
  }  // end of addPostProcessingTask

  void Simulation::addTimeStepValidator(
      const ExternalTimeStepValidator &v) noexcept {
    this->timeStepValidator->addValidator(v);
  }  // end of addTimeStepValidator

  void Simulation::addTimeStepValidator(
      const std::string &n, const ExternalTimeStepValidator &v) noexcept {
    this->timeStepValidator->addValidator(n, v);
  }  // end of addTimeStepValidator

  // bool Simulation::addMonitor(Context &ctx, const
  // std::shared_ptr<AbstractSimulationMonitor> &m) noexcept
  // {
  //   if (m.get() == nullptr) {
  //     return ctx.registerErrorMessage("invalid monitor");
  //   }
  //   this->monitors.push_back(m);
  //   return true;
  // }    // end of addMonitor
  //
  // bool Simulation::addMonitor(Context &ctx, std::string_view n, const
  // Parameters &params) noexcept
  // {
  //   const auto &f = SimulationMonitorFactory::get();
  //   return this->addMonitor(ctx, f.create(ctx, std::string{n},
  //   this->physicalSystem, params));
  // }    // end of addMonitor
  //
  bool Simulation::addTimeIncrementComputer(
      Context &ctx,
      const std::shared_ptr<AbstractTimeIncrementComputer> &m) noexcept {
    if (m.get() == nullptr) {
      return ctx.registerErrorMessage("invalid time increment computer");
    }
    this->timeIncrementComputers.push_back(m);
    return true;
  }  // end of addTimeIncrementComputer

  //   bool Simulation::addTimeIncrementComputer(Context &ctx,
  //                                             std::string_view n,
  //                                             const Parameters &params)
  //                                             noexcept {
  //     const auto &f = TimeIncrementComputerFactory::get();
  //     return this->addTimeIncrementComputer(
  //         ctx, f.create(ctx, std::string{n}, this->physicalSystem, params));
  //   }  // end of addTimeIncrementComputer

  std::pair<ExitStatus, std::optional<SimulationOutput>> Simulation::run(
      Context &ctx) noexcept {
    if (this->timesDescription.size() < 2) {
      std::ignore = ctx.registerErrorMessage(
          "invalid number of times defined for running a simulation");
      return {ExitStatus::unrecoverableError, {}};
    }
    //
    auto s = ExitStatus{};
    auto state = SimulationRunState{};
    //
    auto updateAndSynchronize = [&s](const ExitStatus o) {
      s.update(o);
      s = synchronize(s);
    };
    //
    //   if (!this->physicalSystem.isInitialized()) {
    //     updateAndSynchronize(this->physicalSystem.initialize(ctx));
    //     if (!s.shallContinue()) {
    //       ctx.registerErrorMessage("initialization of the physical system
    //       failed"); return {s};
    //     }
    //   }
    //
    auto p_bts = this->timesDescription.begin();
    auto p_ets = std::next(p_bts);
//    const auto dt = *p_ets - *p_bts;
//   updateAndSynchronize(this->physicalSystem.updateClock(ctx, *p_bts,
//   dt)); if (!s.shallContinue()) {
//     ctx.registerErrorMessage("updating the clock failed");
//     return {s};
//   }
#pragma message("HERE")
    //   for (const auto &[n, t] : this->initializationTasks) {
    //     updateAndSynchronize(invoke(ctx, t));
    //     if (!s.shallContinue()) {
    //       ctx.registerErrorMessage("initialization task '" + n + "'failed");
    //       return {s};
    //     }
    //   }
    //   if (!this->physicalSystem.executeInitialPostProcessingTasks(ctx)) {
    //     return {ExitStatus::unrecoverableError};
    //   }
    // loop over the time sequences
    auto pdt = std::optional<real>{};
    auto pe = this->timesDescription.end();
    auto output = SimulationOutput{};
    //     output["numberOfTimeSteps"] = 0;
    //     output["numberOfSubSteps"] = 0;
    //     output["timeStepOutputs"] = std::vector<Parameters>{};
    while (p_ets != pe) {
      if (this->independentTemporalSequences) {
        pdt = std::optional<real>{};
      }
      this->simulateOverATemporalSequence(ctx, s, output, state, pdt, *p_bts,
                                          *p_ets, std::next(p_ets) == pe);
      updateAndSynchronize(s);
      if (!s.shallContinue()) {
        // set the time at the beginning of the temporal sequence to   the time
        // at the beginning of the time step, so if the run method
        // is called again
        // one starts at the correct time
        //         *p_bts = this->physicalSystem.getClock()
        //                      .getTimeAtTheBeginningOfTheTimeStep();
        return {s, {}};
      }
      p_bts = this->timesDescription.erase(p_bts);
      p_ets = std::next(p_bts);
      pe = this->timesDescription.end();
      if (state.maximumNumberOfTimeStepsReached) {
        // dt is arbitrary, it just have to be a non-zero value
        //    updateAndSynchronize(this->physicalSystem.updateClock(ctx, *p_bts,
        //    dt));
        if (s.shallStop()) {
          return {s, {}};
        } else {
          return {s, output};
        }
      }
    }
    // dt is arbitrary, it just have to be a non-zero value
    //    updateAndSynchronize(this->physicalSystem.updateClock(ctx, *p_bts,
    //    dt));
    if (s.shallStop()) {
      return {s, {}};
    }
    return {s, output};
  }  // end of run

  std::optional<real> Simulation::getEndOfNextTimeStep(
      Context &ctx,
      const real,
      const real se,
      const real t,
      const std::optional<real> pdt) noexcept {
    const auto odt1 = [this, &ctx, se, t]() -> std::optional<real> {
      auto odt = std::optional<real>{};
      for (auto &c : this->timeIncrementComputers) {
        const auto nodt = c->getNextTimeIncrement(ctx, t, se);
        if (isInvalid(nodt)) {
          return {};
        }
        if (isValid(odt)) {
          odt = std::min(*odt, *nodt);
        } else {
          odt = *nodt;
        }
      }
      return odt;
    }();
    if (isInvalid(odt1)) {
      return {};
    }
    auto dt = std::min(*odt1, se - t);
    //      const auto odt2 =
    //      this->physicalSystem.getNextTimeIncrement(ctx, t, se); if
    //      (isInvalid(odt2)) {
    //        return {};
    //      }
    //      auto dt = std::min(std::min(*odt1, *odt2), se - t);
    //
    if (isValid(pdt)) {
      if (this->limitTimeIncrementIncrease) {
        if (dt > *pdt) {
          const auto max_dt =
              (this->maximalTimeIncrementRelativeIncrease) * (*pdt);
          if (dt > max_dt) {
            dt = max_dt;
          }
        }
      }
      if (this->limitTimeIncrementDecrease) {
        if (dt < *pdt) {
          dt = std::max(dt,
                        (this->maximalTimeIncrementRelativeDecrease) * (*pdt));
        }
      }
    }
    //
    if (isValid(this->minimalTimeIncrement)) {
      if (dt < *(this->minimalTimeIncrement)) {
        return ctx.registerErrorMessage(
            "proposed time increment is lower than the minimal one "
            "(" +
            std::to_string(dt) + " < " +
            std::to_string(*(this->minimalTimeIncrement)) + ")");
      }
    }
    if (isValid(this->maximalTimeIncrement)) {
      if (dt > *(this->maximalTimeIncrement)) {
        dt = *(this->maximalTimeIncrement);
      }
    }
    // balance time increment
    if (this->balanceTimeIncrements) {
      const auto rdt = se - t;
      const auto q = std::floor(rdt / (dt));
      const auto r = rdt / (dt)-q;
      if (r < this->minimalRelativeRemainder) {
        // the remainder is very small, so we replace
        // the proposed time increment by rdt / q
        dt = rdt / q;
      } else if (r < this->maximalRelativeRemainder) {
        // we choose to balance the time increment
        dt = rdt / (q + 1);
      } else {
        // the remainder is large enought so that
        // the last time step would be ok
      }
    }
    //
    return t + dt;
  }  // end of getEndOfNextTimeStep

  void Simulation::simulateOverATemporalSequence(
      Context &ctx,
      ExitStatus &s,
      SimulationOutput &output,
      SimulationRunState &state,
      std::optional<real> &pdt,
      const real sb,
      const real se,
      const bool /* lastTemporalSequence */) noexcept {
    // status of the last successful time step
    auto previousStatus = s;
    //
    auto updateAndSynchronize = [&s](const ExitStatus o) {
      s.update(o);
      s = synchronize(s);
    };
    //
    if (ctx.getVerbosityLevel() >= VerbosityLevel::verboseLevel1) {
      ctx.log() << "Computing time sequence from " << sb << " to " << se
                << std::endl;
    }
    // number of substeps
    auto nFailures = size_type{};
    // duration of the temporal sequence
    const auto sdt = se - sb;
    // beginning of the current time step
    auto t = sb;
    // end of the current time step
    auto ote = this->getEndOfNextTimeStep(ctx, sb, se, t, pdt);
    if (isInvalid(ote)) {
      s = ExitStatus::unrecoverableError;
      return;
    }
    //
    auto getTemporalSequenceDescription = [&sb, &se] {
      return "time sequence from " + std::to_string(sb) + " to " +
             std::to_string(se);
    };
    auto getTimeStepDescription = [&getTemporalSequenceDescription, &t, &ote] {
      return "time step from " + std::to_string(t) + " to " +
             std::to_string(*ote) + " in " + getTemporalSequenceDescription();
    };
    //
    auto prepareSubStepping = [this, &ctx, &s, &t, &ote, &updateAndSynchronize,
                               &previousStatus](const real dte) -> bool {
      //
      if (!ote.has_value()) {
        s = ExitStatus::unrecoverableError;
        return ctx.registerErrorMessage("internal error");
      }
      if (!this->allowSubStepping) {
        return false;
      }
      if ((dte < 0) || (std::fpclassify(dte) == FP_ZERO)) {
        s = ExitStatus::recoverableError;
        return ctx.registerErrorMessage("invalid time increment proposal (" +
                                        std::to_string(dte) + ")");
      }
      if (t + dte > *ote) {
        s = ExitStatus::recoverableError;
        return ctx.registerErrorMessage(
            "new time increment proposal is greater than the previous "
            "one (" +
            std::to_string(dte) + " vs " + std::to_string(*ote - t) + ")");
      }
      //
      auto dt = std::max(
          dte, (this->maximalTimeIncrementRelativeDecrease) * (*ote - t));
      //
      if (this->minimalTimeIncrement.has_value()) {
        if (*(this->minimalTimeIncrement) > dt) {
          s = ExitStatus::recoverableError;
          return ctx.registerErrorMessage(
              "new time increment proposal is greater than the minimal "
              "time "
              "increment (" +
              std::to_string(dt) + " vs " +
              std::to_string(*(this->minimalTimeIncrement)) + ")");
        }
      }
      ote = t + dt;
      // restoring the previous status
      s = previousStatus;
      //
      if (isValid(this->nonlinearEvolutionProblem)) {
        this->nonlinearEvolutionProblem->revert();
      }
      //      updateAndSynchronize(this->physicalSystem.revert(ctx));
      return s.shallContinue() ? true : false;
    };
    //
    auto reportMaximumFailureReached = [&ctx, &s, &t,
                                        getTemporalSequenceDescription]() {
      s = ctx.registerErrorMessage(
          "maximum number of failures or time step rejections "
          "reached at "
          "time " +
          std::to_string(t) + " in " + getTemporalSequenceDescription());
    };
    // loop until the end of the temporal sequence
    auto stop = false;
    while (!stop) {
      // saving the previous status
      previousStatus = s;
      if (!((std::abs(sb - t) <
             sdt * 10 * std::numeric_limits<real>::epsilon()) &&
            (std::abs(se - *ote) <
             sdt * 10 * std::numeric_limits<real>::epsilon()))) {
        // the time step does not convert the whole temporal sequence
        if (ctx.getVerbosityLevel() >= VerbosityLevel::verboseLevel1) {
          ctx.log() << "Computing " << getTimeStepDescription() << std::endl;
        }
      }
      updateAndSynchronize(
          this->simulateOverATimeStep(ctx, output, state, t, *ote));
      if (s == ExitStatus::unrecoverableError) {
        return;
      }
      // resolution did not succeed, calling the convergence failure
      // handler
      if (s == ExitStatus::recoverableError) {
#pragma message("HERE")
        //         output["numberOfSubSteps"] =
        //             get<size_type>(output, "numberOfSubSteps") + 1;
        //
        ++nFailures;
        if (nFailures == this->maximumNumberOfFailuresPerTemporalSequence) {
          reportMaximumFailureReached();
          return;
        }
        auto odt =
            this->convergenceFailureHandler->getNewTimeIncrement(ctx, *ote - t);
        if (isInvalid(odt)) {
          return;
        }
        if (!prepareSubStepping(*odt)) {
          return;
        }
        // flushing the error message
        std::ignore = ctx.registerErrorMessage(
            "The coupling scheme did not converge for the " +
            getTimeStepDescription());
        const auto msg = ctx.getErrorMessage();
        if (ctx.getVerbosityLevel() >= VerbosityLevel::verboseLevel1) {
          ctx.log() << msg << std::endl;
        }
        // restoring the previous status
        s = previousStatus;
        continue;
      }
      // validate the time step
      const auto r = this->timeStepValidator->validate(ctx);
      if (isInvalid(r)) {
        s = ctx.registerErrorMessage("The time step validator failed for the " +
                                     getTimeStepDescription());
        return;
      }
      if (!r->isValid) {
        ++nFailures;
        if (nFailures == this->maximumNumberOfFailuresPerTemporalSequence) {
          reportMaximumFailureReached();
          return;
        }
        // flushing the error message
        if (ctx.getVerbosityLevel() >= VerbosityLevel::verboseLevel1) {
          ctx.log() << "The coupling scheme converged for the "
                    << getTimeStepDescription()
                    << ", but this solution was rejected ";
          if (r->reasons.empty()) {
            ctx.log() << "(unspecified reason)";
          } else if (r->reasons.size() == 1) {
            ctx.log() << "for the following reason: " << r->reasons[0];
          } else {
            ctx.log() << "for the following reasons:";
            for (const auto &reason : r->reasons) {
              ctx.log() << "\n- " << reason;
            }
          }
          ctx.log() << std::endl;
        }
        if (!prepareSubStepping(r->recommendedTimeIncrement)) {
          return;
        }
        // restoring the previous status
        s = previousStatus;
        continue;
      }
      // the time step if valid:
      // - here we test if we are at the end of the temporal sequence
      state.timeSinceLastPostProcessing += *ote - t;
      stop =
          std::abs(se - *ote) < sdt * 10 * std::numeric_limits<real>::epsilon();
      auto explicitMarkedPostProcessingTime = stop;
      if (isValid(this->numberOfTimeStepsBetweenPostProcessings)) {
        explicitMarkedPostProcessingTime =
            explicitMarkedPostProcessingTime ||
            ((state.numberOfTimeSteps %
              (*(this->numberOfTimeStepsBetweenPostProcessings))) == 0);
      }
      if (isValid(this->timeBetweenPostProcessings)) {
        if (state.timeSinceLastPostProcessing >
            *(this->timeBetweenPostProcessings)) {
          explicitMarkedPostProcessingTime = true;
          state.timeSinceLastPostProcessing = real{};
        }
      }
      // executing post-processing tasks
      if (isValid(this->nonlinearEvolutionProblem)) {
        this->nonlinearEvolutionProblem->executePostProcessings(t, *ote - t);
      }
      //       if (!this->physicalSystem.executePostProcessingTasks(
      //               ctx, explicitMarkedPostProcessingTime)) {
      //         s = ctx.registerErrorMessage("The time step validator failed
      //         for the " +
      //                                      getTimeStepDescription());
      //         return;
      //       }
#pragma message("HERE")
//       for (const auto &[n, pt] : this->postProcessingTasks) {
//         updateAndSynchronize(invoke(ctx, pt));
//         if (!s.shallContinue()) {
//           ctx.registerErrorMessage("post-processing task '" + n + "'
//           failed"); return;
//         }
//       }
#pragma message("HERE")
      //      // executing monitors
      //       for (auto &m : this->monitors) {
      //         updateAndSynchronize(m->execute(ctx, output,
      //         lastTemporalSequence)); if (!s.shallContinue()) {
      //           ctx.registerErrorMessage("calling a monitor failed");
      //           s = ExitStatus::unrecoverableError;
      //           return;
      //         }
      //       }
      // updating the physical system
      if (isValid(this->nonlinearEvolutionProblem)) {
        this->nonlinearEvolutionProblem->update();
      }
      //      updateAndSynchronize(this->physicalSystem.update(ctx));
      if (!s.shallContinue()) {
        return;
      }
      if (state.maximumNumberOfTimeStepsReached) {
        return;
      }
      // updating the previous status
      previousStatus = s;
      // updating the current time and compute the next time step, if required
      pdt = *ote - t;
      t = *ote;
      if (!stop) {
        ote = this->getEndOfNextTimeStep(ctx, sb, se, t, pdt);
        if (isInvalid(ote)) {
          s = ExitStatus::unrecoverableError;
          return;
        }
      }
    }
  }  // end of simulateOverATemporalSequence

  ExitStatus Simulation::simulateOverATimeStep(Context &ctx,
                                               SimulationOutput &output,
                                               SimulationRunState &state,
                                               const real tb,
                                               const real te) noexcept {
    ++(state.numberOfTimeSteps);
    if (isValid(this->maximumNumberOfTimeSteps)) {
      state.maximumNumberOfTimeStepsReached =
          (state.numberOfTimeSteps == *(this->maximumNumberOfTimeSteps));
    }
    //
    auto s = ExitStatus{};
    //
    //     auto updateAndSynchronize = [&s](const ExitStatus o) {
    //       s.update(o);
    //       s = synchronize(s);
    //     };
    if (ctx.getVerbosityLevel() >= VerbosityLevel::verboseLevel0) {
      ctx.log() << "Computing next state from time " << tb << " to " << te
                << std::endl;
    }
    if (isValid(this->nonlinearEvolutionProblem)) {
      this->nonlinearEvolutionProblem->setup(tb, te - tb);
    }
    //    updateAndSynchronize(this->physicalSystem.updateClock(ctx, tb, te -
    //    tb));
    if (!s.shallContinue()) {
      std::ignore = ctx.registerErrorMessage("updating the clock failed");
      return s;
    }
    //     updateAndSynchronize(
    //         this->physicalSystem.updateLoadingsAtTheBeginningOfTheTimeStep(ctx));
    //     if (!s.shallContinue()) {
    //       ctx.registerErrorMessage(
    //           "updating loading at the beginning of the time step failed");
    //       return s;
    //      }
    //      updateAndSynchronize(this->physicalSystem.performInitializationTaksAtTheBeginningOfTheTimeStep(ctx));
    //     if (!s.shallContinue()) {
    //       ctx.registerErrorMessage(
    //           "performing initialization tasks at the beginning of the time
    //           step " "failed");
    //       return s;
    //     }
    // computing the next time step
    const auto r = [this, &ctx, &tb, &te] {
      if (isValid(this->nonlinearEvolutionProblem)) {
        const auto routput =
            this->nonlinearEvolutionProblem->solve(ctx, tb, te - tb);
        if (isInvalid(routput)) {
          return ExitStatus::recoverableError;
        }
      }
      return ExitStatus::success;
    }();
    // this->physicalSystem.computeNextState(ctx);
    if (r.shallContinue()) {
      for (auto &c : this->timeIncrementComputers) {
        if (!c->prepareNextTimeStep(ctx)) {
          return ExitStatus::unrecoverableError;
        }
      }
      //       if (this->keepOutputs) {
      //         auto timeStepOutput = Parameters{};
      //         timeStepOutput["beginning of time step"] = tb;
      //         timeStepOutput["end of time step"] = te;
      //         timeStepOutput["computeNextStateOutput"] =
      //             static_cast<const Parameters &>(*r);
      //         output["numberOfTimeSteps"] =
      //             get<size_type>(output, "numberOfTimeSteps") + 1;
      //         auto ts = get<std::vector<Parameters>>(output,
      //         "timeStepOutputs"); ts.push_back(timeStepOutput);
      //         output["timeStepOutputs"] = ts;
      //       }
    }
    return r;
  }  // end of simulateOverATimeStep

  bool Simulation::setMaximumNumberOfTimeSteps(Context &ctx,
                                               const size_type n) noexcept {
    if (n < 1) {
      return ctx.registerErrorMessage("invalid maximum number of time step (" +
                                      std::to_string(n) + ")");
    }
    this->maximumNumberOfTimeSteps = n;
    return true;
  }  // end of setMaximumNumberOfTimeSteps

  void Simulation::unsetMaximumNumberOfTimeSteps() noexcept {
    this->maximumNumberOfTimeSteps.reset();
  }

  bool Simulation::setNumberOfTimeStepsBetweenPostProcessings(
      Context &ctx, const size_type n) noexcept {
    if (n < 1) {
      return ctx.registerErrorMessage("invalid post-processing periodicity (" +
                                      std::to_string(n) + ")");
    }
    this->numberOfTimeStepsBetweenPostProcessings = n;
    return true;
  }  // end of setNumberOfTimeStepsBetweenPostProcessings

  void Simulation::unsetNumberOfTimeStepsBetweenPostProcessings() noexcept {
    this->numberOfTimeStepsBetweenPostProcessings.reset();
  }

  bool Simulation::setTimeBetweenPostProcessings(Context &ctx,
                                                 const real dt) noexcept {
    if (dt < std::numeric_limits<real>::min()) {
      return ctx.registerErrorMessage("invalid time between post-processing (" +
                                      std::to_string(dt) + ")");
    }
    this->timeBetweenPostProcessings = dt;
    return true;
  }  // end of setTimeBetweenPostProcessings

  void Simulation::unsetTimeBetweenPostProcessings() noexcept {
    this->timeBetweenPostProcessings.reset();
  }

  const Simulation::TimesDescription &Simulation::getTimes() const noexcept {
    return this->timesDescription;
  }

  // static std::pair<size_type, const ResourcesUsageReport *>
  // findSimulationReport(const size_type depth, const
  // std::vector<ResourcesUsageReport> &profilings) noexcept
  // {
  //   for (const auto &p : profilings) {
  //     if (p.getLabel() == "simulation::run") {
  //       return {depth, &p};
  //     }
  //   }
  //   auto mindepth   = depth;
  //   auto *profiling = static_cast<const ResourcesUsageReport *>(nullptr);
  //   for (const auto &p : profilings) {
  //     const auto [d, ptr] = findSimulationReport(depth + 1, p.profilings);
  //     if (ptr != nullptr) {
  //       if (profiling == nullptr) {
  //         mindepth  = d;
  //         profiling = ptr;
  //       }
  //       else if (d < mindepth) {
  //         mindepth  = d;
  //         profiling = ptr;
  //       }
  //     }
  //   }
  //   return {mindepth, profiling};
  // }    // end of findSimulationReport
  //
  // std::optional<const ResourcesUsageReport *> findSimulationReport(const
  // std::vector<ResourcesUsageReport> &profilings) noexcept
  // {
  //   const auto r = findSimulationReport(0, profilings);
  //   if (r.second == nullptr) {
  //     return {};
  //   }
  //   return r.second;
  // }    // end of findSimulationReport

}  // end of namespace mfem_mgis
