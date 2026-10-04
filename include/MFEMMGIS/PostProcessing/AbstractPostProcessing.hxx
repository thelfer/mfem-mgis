/*!
 * \file   MFEMMGIS/PostProcessing/AbstractPostProcessing.hxx
 * \brief  This file declares the AbstractPostProcessing class
 * \author Thomas Helfer
 * \date   27/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTPOSTPROCESSING_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTPOSTPROCESSING_HXX

#include <string>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStep.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct PhysicalSystem;
  struct Context;

  //! \brief an abstract class describing a post-processing
  struct MFEM_MGIS_EXPORT AbstractPostProcessing {
    //! \return the name of the post-processing
    [[nodiscard]] virtual std::string getName() const noexcept = 0;
    //! \return the underlying physical system
    [[nodiscard]] virtual PhysicalSystem &getPhysicalSystem() noexcept = 0;
    //! \return the underlying physical system
    [[nodiscard]] virtual const PhysicalSystem &getPhysicalSystem()
        const noexcept = 0;
    //! \return the `executeInitialPostProcessing` has already been called
    [[nodiscard]] virtual bool
    hasExecuteInitialPostProcessingTasksAlreadyBeenCalled() const noexcept = 0;
    /*!
     * \brief execute the post-processing at the beginning of the simulation.
     * For instance, this method may display the initial values of the state
     * variables.
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] virtual bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept = 0;
    /*!
     * \brief execute the post-processing at the end of a time step, after
     * convergence.
     *
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] isPostProcessingRequired: boolean stating if the time at
     * the end of the time step is a post-processing time. By default, the end
     * of a temporal sequence is a post-processing time (but the definition of
     * other post-processing times is possible). The implementations shall let
     * the user choose if the post-processing must be executed at every time
     * step or only at post-processing times by accepting an `AllTimeSteps`
     * parameter. The default behavior is left to the implementation, but the
     * following rule of thumb is that a lightweight post-processing shall be
     * executed by default at each end of time steps and that a heavy
     * post-processing (in execution time and/or size of the generated file(s))
     * shall be executed at post-processing times by default.
     * \return true on success
     */
    [[nodiscard]] virtual bool executePostProcessingTasks(
        Context &ctx,
        const TimeStep &ts,
        const bool isPostProcessingRequired) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractPostProcessing() noexcept;
  };  // end of  class AbstractPostProcessing

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTPOSTPROCESSING_HXX */