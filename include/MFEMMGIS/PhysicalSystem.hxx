/*!
 * \file   MFEMMGIS/PhysicalSystem.hxx
 * \brief  This file declares the `PhysicalSystem` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_PHYSICAL_SYSTEM_HXX
#define LIB_MFEM_MGIS_PHYSICAL_SYSTEM_HXX

#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ExitStatus.hxx"
#include "MFEMMGIS/ComputeNextStateOutput.hxx"

namespace mfem_mgis {

  //! forward declarations
  struct TimeStep;
  struct Parameters;
  struct AbstractCouplingScheme;
  struct AbstractModel;
  struct AbstractPostProcessing;

  /*!
   * \brief the physical system class is meant to gather all the
   * information required to describe the evolution of a physical
   * system over a time step.
   */
  struct MFEM_MGIS_EXPORT PhysicalSystem {
    /*!
     * \brief constructor
     */
    PhysicalSystem();
    /*!
     * \return a description of the physical system
     * \param[in] ctx: execution context
     * \param[in] b: boolean being the default value for information requests.
     * \param[in] parameters: dictionary that allows the parametrize the output.
     *
     * This dictionary is meant to contains boolean values corresponding
     * to information requests. The default value of the boolean is given by the
     * `b` parameter.
     */
    std::optional<std::string> describe(Context &,
                                        const bool,
                                        const Parameters &) const noexcept;
    //! \return if the coupling scheme is defined
    bool isCouplingSchemeDefined() const noexcept;
    /*!
     * \brief set the coupling scheme
     * \param[out] ctx: execution context.
     * \param[in]  c: coupling scheme.
     */
    [[nodiscard]] bool setCouplingScheme(
        Context &, std::shared_ptr<AbstractCouplingScheme>) noexcept;
    //     /*!
    //      * \brief set the coupling scheme
    //      * \param[out] ctx: execution contex
    //      * \param[in]  n: name of the model
    //      * \param[in]  p: parameters used to initialize the coupling scheme
    //      */
    //     [[nodiscard]] bool setCouplingScheme(Context &,
    //                                          std::string_view,
    //                                          const Parameters &) noexcept;
    /*!
     * \brief set the unique model. This methods sets a default
     * coupling scheme that only call this model once (per time step).
     *
     * \param[out] ctx: execution context.
     * \param[in]  m: model.
     */
    [[nodiscard]] bool setModel(Context &,
                                std::shared_ptr<AbstractModel>) noexcept;
    //     /*!
    //      * \brief set the unique model. This methods sets a default
    //      * coupling scheme that only call this model once (per time step).
    //      *
    //      * \param[out] ctx: execution context
    //      * \param[in]  n: name of the model
    //      * \param[in]  p: parameters used to initialize the model
    //      */
    //     [[nodiscard]] bool setModel(Context &,
    //                                 std::string_view,
    //                                 const Parameters &) noexcept;
    /*!
     * \brief add a new post-processing
     * \param[in] ctx: execution context
     * \param[in] n: name of the post-processing
     * \param[in] parameters: parameters passed to the post-processing
     */
    [[nodiscard]] bool addPostProcessing(Context &,
                                         std::string_view,
                                         const Parameters &) noexcept;
    /*!
     * \brief add a new post-processing
     * \param[in] ctx: execution context
     * \param[in] p: post-processing
     */
    [[nodiscard]] bool addPostProcessing(
        Context &, std::shared_ptr<AbstractPostProcessing>) noexcept;
    /*!
     * \brief update loadings.
     *
     * \param[in, out] ctx: execution contex
     * \param[in] ts: description of the time step
     *
     * \note This method must be called at the beginning of the time step.
     */
    [[nodiscard]] bool updateLoadingsAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept;
    /*!
     * \brief perform initialization tasks.
     *
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     *
     * \note this method shall be called after
     * `updateLoadingsAtTheBeginningOfTheTimeStep`
     * \note this method shall be called before `computeNextState`
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept;
    /*!
     * \brief This method is called at the beginning of a time step to determine
     * a suitable time increment.
     *
     * \return the next time increment.
     * \param[in] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    [[nodiscard]] std::optional<real> getNextTimeIncrement(
        Context &, const real, const real) const noexcept;
    /*!
     * \brief compute the next state
     * \param[out] ctx: execution context
     * \param[in] ts: description of the time step
     *
     * \note this method shall be called after `computeNextState`
     * \note this method shall be called before `update`
     */
    std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept;
    /*!
     * \brief execute post-processings at the beginning of the simulation. For
     * instance, this method may display the initial values of the state
     * variables.
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time step
     */
    [[nodiscard]] bool executeInitialPostProcessingTasks(Context &,
                                                         const real) noexcept;
    /*!
     * \brief execute post-processings at the end of a time step, after
     * convergence.
     *
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] b: boolean stating that if the time at the end of the time
     * step is a post-processing time.
     */
    [[nodiscard]] bool executePostProcessingTasks(Context &,
                                                  const TimeStep &,
                                                  const bool);
    /*!
     * \brief update the system
     * This method shall be called after `computeNextState`
     */
    [[nodiscard]] bool update(Context &) noexcept;
    /*!
     * \brief revert the system to its state at the beginning of the time step
     */
    [[nodiscard]] bool revert(Context &) noexcept;
    //! \brief destructor
    ~PhysicalSystem() noexcept;

   private:
    //! \brief coupling scheme
    std::shared_ptr<AbstractCouplingScheme> coupling_scheme;
    //! \brief list of registered post processing
    std::vector<std::shared_ptr<AbstractPostProcessing>> post_processings;
  };  // end of class PhysicalSystem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PHYSICAL_SYSTEM_HXX */