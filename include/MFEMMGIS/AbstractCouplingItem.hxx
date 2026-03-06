/*!
 * \file   MFEMMGIS/AbstractCouplingItem.hxx
 * \brief  This file declares the `AbstractCouplingItem` class
 * \date   06/12/2022
 */

#ifndef LIB_MFEM_MGIS_ABSTRACT_COUPLING_ITEM_HXX
#define LIB_MFEM_MGIS_ABSTRACT_COUPLING_ITEM_HXX

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <ostream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ExitStatus.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/ComputeNextStateOutput.hxx"

namespace mfem_mgis {

  // forward declaration
  struct TimeStep;
  struct Parameters;

  /*!
   * \brief interface of all coupling items, i.e. objects that can be added to a
   * coupling scheme
   *
   * \note the `getVerbosityLevel` and `getLogStreamPointer` methods
   * are meant to be used by coupling schemes to modify the context passed to
   * the coupling item automatically. There is probably no reason to
   * for other objects to call those methods.
   */
  struct MFEM_MGIS_EXPORT AbstractCouplingItem {
    //! \brief name of the parameter associated with the verbosity level
    static const std::string verbosityLevelParameter;
    //! \brief name of the parameter associated with an output file
    static const std::string logFileParameter;
    //! \return a name describing the coupling item
    [[nodiscard]] virtual std::string getName() const noexcept = 0;
    //! \brief return the mesh discretization
    virtual MeshDiscretization getMeshDiscretization() const noexcept = 0;
    //! \return the list of locations on which the item works
    [[nodiscard]] virtual std::vector<std::string> getLocations()
        const noexcept = 0;
    //! \return the verbosity level to be used when the coupling item is called
    [[nodiscard]] virtual VerbosityLevel getVerbosityLevel() const noexcept = 0;
    /*!
     * \brief associate a log stream to the current item. The given pointer can
     * be null. \param[in] s: log stream
     */
    virtual void setLogStream(std::shared_ptr<std::ostream>) noexcept = 0;
    /*!
     * \return a pointer to a log stream. This pointer may be null if the
     * coupling item does not declare a specific log stream. \param[in] s: log
     * stream
     */
    [[nodiscard]] virtual std::shared_ptr<std::ostream>
    getLogStreamPointer() noexcept = 0;
    /*!
     * \brief set the verbosity level to be used when the coupling item is
     * called \param[in] l: the new verbose level
     */
    virtual void setVerbosityLevel(const VerbosityLevel) noexcept = 0;
    /*!
     * \return a description of the coupling item
     * \param[in] b: boolean being the default value for information requests.
     * \param[in] parameters: dictionary that allows the parametrize the output.
     * This dictionary is meant to contains boolean values corresponding
     * to information requests. The default value of the boolean is given by the
     * `b` parameter, except for the `shortDescription` request which default
     * value is true.
     *
     * \note the short description is excepted to stand in a single line.
     */
    [[nodiscard]] virtual std::optional<std::string> describe(
        Context &, const bool, const Parameters &) const noexcept = 0;
    //     /*!
    //      * \brief declare the minimal depencies of the item.
    //      * \param[in] ctx: execution context
    //      * \param[in] dm: dependencies manager
    //      */
    //     virtual bool declareDependencies(Context &,
    //                                      DependenciesManager &) const
    //                                      noexcept = 0;
    //     /*!
    //      * \brief perform initalize tasks before the allocation of the
    //      resources
    //      * (indeed this method shall declare all the required resources
    //      * (formulations, temporary fields) required to perform the
    //      computations).
    //      * \param[in] ctx: execution context
    //      * \param[in] vf: factory of value evaluators
    //      * \param[in] nf: factory of nodal evaluators
    //      * \param[in] f: factory of evaluators at integration points
    //      * \note this method must be called once every dependencies has been
    //      * resolved. \note this method must be called **before** the
    //      initialization
    //      * of the resources manager of the physical system.
    //      */
    //     virtual bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept = 0;
    /*!
     * \brief method called at the beginning of a time step
     *
     * \param[in,out] ctx: execution context
     * \param[in] ts: description of the time step
     *
     * \note if required, the time step can be retrieved from the clock hold by
     * the physical system
     */
    [[nodiscard]] virtual bool
    performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept = 0;
    /*!
     * \brief This method is called at the beginning of a time step to determine
     * a suitable time increment.
     *
     * \return the next time increment.
     * \param[in, out] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    [[nodiscard]] virtual std::optional<real> getNextTimeIncrement(
        Context &, const real, const real) const noexcept = 0;
    /*!
     * \brief compute the state of the system at the end of of the time step
     *
     * \param[in,out] ctx: execution context
     * \param[in] ts: description of the time step
     *
     * The returned structure contains the information about the execution.
     * Each item must document what is returned. Developpers are advised to
     * return uniform information for items having the same role (i.e. coupling,
     * models, etc...)
     *
     * For couplings, we expect the following information:
     *
     * - `NumberOfIterations` (integer): the number of iterations up to
     * convergence
     * - `IterationsOutputs` (Vector of parameters): by iteration, the outputs
     *   of each coupling item. Each element of the vector corresponds to an
     *   iteration, the first element corresponding to the first iteration. Each
     *   element of this vector is a vector of dictionaries whose structure is
     *   given below.
     * - `ItemsOutputs` (vector of dictionaries): the outputs of each coupling
     * item at the last iteration.
     *
     * For an iteration, the outputs of each coupling item of this vector must
     * contains
     *   - `Name` (string): the name of item
     *   - `Description` (string): the description of the item
     *   - `Output` (dictionary): the value returned by `computeNextState` for
     * the last iteration
     *
     * For models using an implicit solver, we expect the following output:
     *
     * - `Solver` (dictionary): this dictionary may contain the following
     * entries:
     *   - `NumberOfIterations` (integer): the number of iterations to reach
     * convergence
     *
     * \note if required, the time step can be retrieved from the clock hold by
     * the physical system
     */
    [[nodiscard]] virtual std::pair<ExitStatus,
                                    std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept = 0;
    /*!
     * \brief execute post-processings at the beginning of the simulation. For
     * instance, this method may display the initial values of the state
     * variables.
     *
     * \param[in,out] ctx: execution context
     * \param[in] t: initial time
     *
     * \note if required, the time at the beginning of the time step can be
     * retrieved from the clock hold by the physical system
     */
    [[nodiscard]] virtual bool executeInitialPostProcessingTasks(
        Context &, const real) noexcept = 0;
    /*!
     * \brief execute post-processings at the end of a time step, after
     * convergence.
     *
     * \param[in,out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] b: boolean stating that if the time at the end of the time
     * step is a post-processing time.
     *
     * \note if required, the time at the end of the time step can be retrieved
     * from the clock hold by the physical system
     */
    [[nodiscard]] virtual bool executePostProcessingTasks(
        Context &, const TimeStep &, const bool) noexcept = 0;
    /*!
     * \brief update the state of the system for the next time step
     *
     * This method is called at the end of the time step once the state of
     * the system at the end of the time step is known. This state is generally
     * defined by a set of unknown fields (usually defined as node fields or
     * element fields) and a set of internal state variables (usually defined as
     * fields at integration points).
     *
     * This method is first meant to copy the fields at the end of the time on
     * the fields at the beginning of the time step.
     */
    [[nodiscard]] virtual bool update(Context &) noexcept = 0;
    /*!
     * \brief revert the state of the system at the beginning of the time step
     *
     * This method is typically called in case of non convergence of the
     * coupling scheme to restart the computation with a smaller the time step.
     *
     * This method is first meant to copy the fields at the end of the time on
     * the fields at the beginning of the time step.
     */
    [[nodiscard]] virtual bool revert(Context &) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractCouplingItem();
  };

  /*!
   * \return a description of the given item from information returned by the
   * `getName` and `getMeshSetsNames` methods \param[in] i: coupling item
   */
  [[nodiscard]] MFEM_MGIS_EXPORT std::string getShortDescription(
      const AbstractCouplingItem &) noexcept;
  /*!
   * \return a description of the parameters that shall be used by all coupling
   * items. Those parameters are related to the customization of a `Context`
   * (verbosity level and/or log stream).
   */
  [[nodiscard]] MFEM_MGIS_EXPORT std::map<std::string, std::string>
  getCouplingItemParametersDescription() noexcept;
  /*!
   * \brief handle parameters common to all coupling items.
   * Those parameters are related to the customization of a `Context` (verbosity
   * level and/or log stream). \param[in] ctx: execution context \param[in] i:
   * coupling item \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT bool handleCouplingItemParameters(Context &,
                                                     AbstractCouplingItem &,
                                                     const Parameters &);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACT_COUPLING_ITEM_HXX */
