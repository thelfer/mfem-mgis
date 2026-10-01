/*!
 * \file   MFEMMGIS/Simulation.hxx
 * \brief  This file declares the Simulation class.
 * \date   15/05/2023
 */

#ifndef LIB_MFEM_MGIS_SIMULATION_HXX
#define LIB_MFEM_MGIS_SIMULATION_HXX

#include <map>
#include <span>
#include <vector>
#include <algorithm>
#include <functional>
#include <initializer_list>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/ExitStatus.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declarations
  struct PhysicalSystem;
  struct AbstractNonLinearEvolutionProblem;
  struct AbstractSimulationMonitor;
  struct AbstractConvergenceFailureHandler;
  struct AbstractTimeIncrementComputer;
  struct AbstractTimeStepValidator;

  /*!
   * \brief a structure returned by the run method of the Simulation struct.
   *
   * This structure is mostly a wrapper around the `Parameters`
   * struct.
   */
  struct MFEM_MGIS_EXPORT SimulationOutput : public Parameters {
    using Parameters::Parameters;
    using Parameters::operator=;
    //! \brief default constructor
    SimulationOutput() noexcept;
    //! \brief move constructor
    SimulationOutput(SimulationOutput &&) noexcept;
    //! \brief copy constructor
    SimulationOutput(const SimulationOutput &) noexcept;
    //! \brief move assignment
    SimulationOutput &operator=(SimulationOutput &&) noexcept;
    //! \brief copy assignment
    SimulationOutput &operator=(const SimulationOutput &) noexcept;
    //! \brief destructor
    ~SimulationOutput() noexcept;
  };  // end of SimulationOutput

  /*!
   * \brief the simulation struct handles the computation
   * of the evolution of the state of a physical system
   * over a given number of time sequences.
   */
  struct MFEM_MGIS_EXPORT Simulation {
    /*!
     * \brief a helper struct to describe a set of time steps
     */
    struct MFEM_MGIS_EXPORT TimesDescription : private std::vector<real> {
      /*!
       * \brief constructor
       * \param[in] tb: time at the beginning of the simulation.
       * \param[in] te: final time.
       */
      TimesDescription(const real tb, const real te);
      /*!
       * \brief constructor
       * \param[in] tb: time at the beginning of the simulation.
       * \param[in] te: final time.
       * \param[in] n: number of temporal sequences.
       */
      TimesDescription(const real tb, const real te, const size_type n);
      /*!
       * \brief constructor
       * \param[in] times: list of times. The times must be ordered.
       */
      TimesDescription(std::span<const real> times);
      //! \brief move constructor
      TimesDescription(TimesDescription &&) noexcept;
      //! \brief copy constructor
      TimesDescription(const TimesDescription &) noexcept;
      //! \brief move assignment
      TimesDescription &operator=(TimesDescription &&) noexcept;
      //! \brief copy assignment
      TimesDescription &operator=(const TimesDescription &) noexcept;
      //
      using std::vector<real>::const_iterator;
      using std::vector<real>::cbegin;
      using std::vector<real>::cend;
      using std::vector<real>::size;
      //! \brief destructor
      ~TimesDescription() noexcept;

     private:
      // this allows the Simulation struct to access the erase method
      friend struct Simulation;
    };
    /*!
     * \brief external time step validator. It returns a boolean stating if
     * the time step is valid and a recommended time increment.
     */
    using ExternalTimeStepValidator = std::function<std::pair<bool, real>()>;
    //! \brief task executed at the beginning of the simulation
    using InitializationTask = std::function<void()>;
    //! \brief task executed at the end of each successful time step
    using PostProcessingTask = std::function<void()>;
    //! \return a description of the parameters of a simulation
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] parameters: parameters passed to the simulation
     *
     * The main parameters are:
     *
     * - `Times`: list of times defining temporal sequences.
     * - `TimeIncrementComputer`: strategy used to determine the next time
     *   step.
     * - `ConvergenceFailureHandler`: strategy used to determine how a
     *   convergence failure of the resolution shall be handled.
     * - `TimeStepValidator`: strategy used to validate the time step.
     */
    Simulation(Context &ctx,
               AbstractNonLinearEvolutionProblem &p,
               const Parameters &parameters);

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] times: description of the temporal sequences
     */
    Simulation(Context &ctx,
               AbstractNonLinearEvolutionProblem &p,
               const TimesDescription &times);

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] times: list of times
     */
    Simulation(Context &ctx,
               AbstractNonLinearEvolutionProblem &p,
               const std::initializer_list<real> &times);

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] parameters: parameters passed to the simulation
     *
     * The main parameters are:
     * - `Times`: list of times defining temporal sequences.
     * - `TimeIncrementComputer`: strategy used to determine the next
     *   time step.
     * - `ConvergenceFailureHandler`: strategy used to determine how a
     *   convergence failure of the resolution shall be handled.
     * - `TimeStepValidator`: strategy used to validate the time step.
     */
    Simulation(Context &ctx, PhysicalSystem &ps, const Parameters &parameters);

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] times: description of the temporal sequences
     */
    Simulation(Context &ctx, PhysicalSystem &ps, const TimesDescription &times);

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] times: list of times
     */
    Simulation(Context &ctx,
               PhysicalSystem &ps,
               const std::initializer_list<real> &times);
    /*!
     * \brief add a task to be performed at the beginning of the simulation
     * \param[in] n: name of the initialization task
     * \param[in] t: initialization task
     */
    void addInitializationTask(const std::string &n,
                               const InitializationTask &t) noexcept;
    /*!
     * \brief add a task to be performed at the beginning of the simulation
     * \param[in] t: initialization task
     */
    void addInitializationTask(const InitializationTask &t) noexcept;
    /*!
     * \brief add a task to be performed at the end of a time step
     * \param[in] n: name of the post-processing task
     * \param[in] t: post-processing task
     */
    void addPostProcessingTask(const std::string &n,
                               const PostProcessingTask &t) noexcept;
    /*!
     * \brief add a task to be performed at the end of a time step
     * \param[in] t: post-processing task
     */
    void addPostProcessingTask(const PostProcessingTask &t) noexcept;
    /*!
     * \brief add a user defined time step validator
     * \param[in] v: external validator
     */
    void addTimeStepValidator(const ExternalTimeStepValidator &v) noexcept;
    /*!
     * \brief add a user defined time step validator
     * \param[in] n: name of the external validator
     * \param[in] v: external validator
     */
    void addTimeStepValidator(const std::string &n,
                              const ExternalTimeStepValidator &v) noexcept;
    /*!
     * \brief add a new time increment computer
     * \param[in, out] ctx: execution context
     * \param[in] m: time increment computer
     * \return true on success
     */
    [[nodiscard]] bool addTimeIncrementComputer(
        Context &ctx,
        const std::shared_ptr<AbstractTimeIncrementComputer> &m) noexcept;
    //   /*
    //    * \brief add a new time increment computer
    //    * \param[in, out] ctx: execution context
    //    * \param[in] n: name
    //    * \param[in] params: parameters passed to the time increment
    //    * computer
    //    * \return true on success
    //    */
    //   [[nodiscard]] bool addTimeIncrementComputer(Context &,
    //   std::string_view, const Parameters &) noexcept;
    /*!
     * \brief add a new monitor
     * \param[in, out] ctx: execution context
     * \param[in] m: monitor
     * \return true on success
     */
    [[nodiscard]] bool addMonitor(
        Context &ctx,
        const std::shared_ptr<AbstractSimulationMonitor> &m) noexcept;
    /*!
     * \brief add a new monitor
     * \param[in, out] ctx: execution context
     * \param[in] n: name
     * \param[in] params: parameters passed to the monitor
     * \return false, this method is not implemented yet
     */
    [[nodiscard]] bool addMonitor(Context &ctx,
                                  std::string_view n,
                                  const Parameters &params) noexcept;
    /*!
     * \brief set the maximum number of time steps
     * \param[in, out] ctx: execution context
     * \param[in] n: maximum number of time steps
     * \return true on success
     */
    [[nodiscard]] bool setMaximumNumberOfTimeSteps(Context &ctx,
                                                   const size_type n) noexcept;
    //! \brief clear the maximum number of time steps
    void unsetMaximumNumberOfTimeSteps() noexcept;
    /*!
     * \brief set the post-processing periodicity
     * \param[in, out] ctx: execution context
     * \param[in] n: number of time steps between two post-processings
     * \return true on success
     */
    [[nodiscard]] bool setNumberOfTimeStepsBetweenPostProcessings(
        Context &ctx, const size_type n) noexcept;
    //! \brief clear the post-processing periodicity
    void unsetNumberOfTimeStepsBetweenPostProcessings() noexcept;
    /*!
     * \brief set the time between post-processings
     * \param[in, out] ctx: execution context
     * \param[in] dt: time between two post-processings
     * \return true on success
     */
    bool setTimeBetweenPostProcessings(Context &ctx, const real dt) noexcept;
    //! \brief clear the time between post-processings
    void unsetTimeBetweenPostProcessings() noexcept;
    //! \return the current time description
    [[nodiscard]] const TimesDescription &getTimes() const noexcept;
    /*!
     * \brief run the simulation
     * \param[in, out] ctx: execution context
     * \return the exit status and, on success, the outputs of the simulation
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<SimulationOutput>> run(
        Context &ctx) noexcept;

   protected:
    /*!
     * \brief initialize the simulation from parameters
     *
     * \param[in, out] ctx: execution context
     * \param[in] parameters: parameters used to initialize the simulation
     */
    void treatParameters(Context &ctx, const Parameters &parameters);
    /*!
     * \brief structure describing the state of a simulation run
     */
    struct SimulationRunState {
      //! \brief number of time steps done
      size_type numberOfTimeSteps = size_type{};
      //! \brief time since the last post-processing
      real timeSinceLastPostProcessing = real{};
      //! \brief boolean stating if the maximum number of time steps is reached
      bool maximumNumberOfTimeStepsReached = false;
    };
    /*!
     * \brief return the end of the next time step
     * \param[in, out] ctx: execution context
     * \param[in] sb: beginning of the time sequence.
     * \param[in] se: end of the time sequence.
     * \param[in] t: current time.
     * \param[in] pdt: previous time increment, if defined.
     * \return the end of the next time step
     */
    [[nodiscard]] std::optional<real> getEndOfNextTimeStep(
        Context &ctx,
        const real sb,
        const real se,
        const real t,
        const std::optional<real> pdt) noexcept;
    /*!
     * \brief run the simulation over a temporal sequence
     * \param[in, out] ctx: execution context
     * \param[in, out] s: current status of the simulation.
     * \param[in, out] output: output of the simulation.
     * \param[in, out] state: simulation state.
     * \param[out] t: time at the end of the last successful increment or time
     * at the beginning of the time step if no successful time step has been
     * made.
     * \param[in, out] pdt: previous time increment.
     * \param[in] sb: beginning of the time sequence.
     * \param[in] se: end of the time sequence.
     * \param[in] lastTemporalSequence: boolean stating if it is the last
     * temporal sequence
     */
    void simulateOverATemporalSequence(
        Context &ctx,
        ExitStatus &s,
        SimulationOutput &output,
        SimulationRunState &state,
        real &t,
        std::optional<real> &pdt,
        const real sb,
        const real se,
        const bool lastTemporalSequence) noexcept;
    /*!
     * \brief run the simulation over a time step
     * \param[in, out] ctx: execution context
     * \param[in, out] output: output of the simulation.
     * \param[in, out] state: simulation state.
     * \param[in] ts: description of the time step
     * \return the exit status of the time step
     */
    ExitStatus simulateOverATimeStep(Context &ctx,
                                     SimulationOutput &output,
                                     SimulationRunState &state,
                                     const TimeStep &ts) noexcept;
    //! \brief the physical system
    OptionalReference<PhysicalSystem> physicalSystem;
    //! \brief evolution problem solved
    OptionalReference<AbstractNonLinearEvolutionProblem>
        nonlinearEvolutionProblem;
    //! \brief the list of temporal sequences
    TimesDescription timesDescription;
    //! \brief initialization tasks
    std::vector<std::pair<std::string, InitializationTask>> initializationTasks;
    //! \brief post-processing tasks
    std::vector<std::pair<std::string, PostProcessingTask>> postProcessingTasks;
    //! \brief monitors
    std::vector<std::shared_ptr<AbstractSimulationMonitor>> monitors;
    /*!
     * \brief boolean stating if sub stepping in case of convergence
     * failure is allowed
     */
    bool allowSubstepping = true;
    /*!
     * \brief boolean stating if the temporal sequences are independent.
     *
     * Currently, this boolean only chooses if the last time increment of the
     * previous temporal sequence is taken into account to bound the first
     * time increment of the current temporal sequence
     * (bounding occurs if this boolean is false)
     */
    bool independentTemporalSequences = true;
    //! \brief maximum number of failures per temporal sequence
    size_type maximumNumberOfFailuresPerTemporalSequence = 10;
    //! \brief minimal time increment
    std::optional<real> minimalTimeIncrement;
    //! \brief maximal time increment
    std::optional<real> maximalTimeIncrement;
    /*!
     * \brief boolean stating if the next time increment is bounded by the
     * previous time increment multiplied by the
     * `maximalTimeIncrementRelativeIncrease` coefficient
     */
    bool limitTimeIncrementIncrease = true;
    /*!
     * \brief boolean stating if the next time increment is bounded below by
     * the previous time increment multiplied by the
     * `maximalTimeIncrementRelativeDecrease` coefficient
     */
    bool limitTimeIncrementDecrease = true;
    /*!
     * \brief coefficient used to determine the maximum ratio between the
     * next time step and the previous one
     */
    real maximalTimeIncrementRelativeIncrease = 1.1;
    /*!
     * \brief coefficient used to determine the minimal ratio between the
     *    next time step and the previous one
     */
    real maximalTimeIncrementRelativeDecrease = 0.2;
    /*!
     * \brief boolean stating if the Simulation struct shall try to
     * balance the time increments
     */
    bool balanceTimeIncrements = true;
    /*!
     * \brief minimal relative remainder for the default time increment
     * balancer
     */
    real minimalRelativeRemainder = 5e-2;
    /*!
     * \brief maximal relative remainder for the default time increment
     * balancer
     */
    real maximalRelativeRemainder = 0.8;
    /*!
     * \brief object called to determine what time step shall be used
     * to restart the computation in case of convergence failure.
     *
     * \note the default handler divides the time increment by two
     */
    std::shared_ptr<AbstractConvergenceFailureHandler>
        convergenceFailureHandler;
    /*!
     * \brief objects called at the beginning of the time step to determine
     * the next time step in the current sequence.
     */
    std::vector<std::shared_ptr<AbstractTimeIncrementComputer>>
        timeIncrementComputers;
    /*!
     * \brief object called to validate the time step
     * after the convergence of the coupling scheme.
     *
     * \note the default validator only calls the external validators.
     */
    std::shared_ptr<AbstractTimeStepValidator> timeStepValidator;
    /*!
     * \brief maximum number of time steps per call to `run`
     *
     * Once this maximum number of time steps is reached, the simulation
     * is stopped without error.
     */
    std::optional<size_type> maximumNumberOfTimeSteps;
    /*!
     * \brief the number of time steps between two post-processings marked
     *  as 'explicitly requested by the user'.
     *
     * Every time a post-processing is called, it receives a boolean
     * stating if the post-processing time has been
     * explicitly requested by the user. Lightweight post-processings
     * may ignore this boolean. This boolean is mostly important for heavy
     *  post-processings in terms of memory, disk usage
     * or computations, such as the `ParaviewExport` post-processing.
     *
     * \note each end of a temporal sequence is always flagged as
     * 'explicitly requested by the user', independently of the
     * `numberOfTimeStepsBetweenPostProcessings` parameter.
     * \note the number of computed time steps is reset at each call of the
     * `run` method. The parameter `numberOfTimeStepsBetweenPostProcessings`
     * thus has no effect if it is greater than `maximumNumberOfTimeSteps`
     * parameter.
     */
    std::optional<size_type> numberOfTimeStepsBetweenPostProcessings;
    /*!
     * \brief the time between two post-processings marked as 'explicitly
     *     requested by the user'.
     *
     * Every time a post-processing is called, it receives a boolean
     * stating if the post-processing time has been
     * explicitly requested by the user. Lightweight post-processings
     * may ignore this boolean. This boolean is mostly important for heavy
     * post-processings in terms of memory, disk usage
     * or computations, such as the `ParaviewExport` post-processing.
     *
     * \note each end of a temporal sequence is always flagged as
     * 'explicitly requested by the user', independently of the
     * `timeBetweenPostProcessings` parameter.
     * \note the time since the last post-processing is reset at each call of
     * the `run` method.
     */
    std::optional<real> timeBetweenPostProcessings;
    /*!
     * \brief keep the outputs of all time steps
     */
    const bool keepOutputs;

   private:
    /*!
     * \brief method called at the end of constructors to define
     * a default time step validator, a default time step computer
     * or a default convergence failure handler, if they were not
     * already defined.
     */
    void completeInitialization();
  };

  // /*
  //  * \brief find in the given profiling reports the first one with the
  //  * label "simulation::run", if any
  //  * \param[in] profilings: profiling reports
  //  * \return the report found, if any
  //  */
  // MFEMMGIS_EXPORT std::optional<const ResourcesUsageReport *>
  // findSimulationReport(
  //     const std::vector<ResourcesUsageReport> &) noexcept;

}  // end of namespace mfem_mgis

#endif LIB_MFEM_MGIS_SIMULATION_HXX
