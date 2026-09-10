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
    /*!
     * \brief execute the post-processing at the beginning of the simulation.
     * For instance, this method may display the initial values of the state
     * variables.
     *
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time step
     *
     * \note if required, the time at the beginning of the time step can be
     * retrieved from the clock hold by the physical system
     */
    [[nodiscard]] virtual bool executeInitialPostProcessingTasks(
        Context &, const real) noexcept = 0;
    /*!
     * \brief execute the post-processing at the end of a time step, after
     * convergence.
     *
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] isPostProcessingRequired: boolean stating that if the time at
     * the end of the time step is a post-processing time. By default, the end
     * of a temporal sequence is a post-processing time (but the definition of
     * other post-processing times is possible). The implementations shall let
     * the user choose if the post-processing must executed at every time steps
     * or only at post-processing times by accepting an `allTimes` parameter.
     * The default behavior is left to the implementation, but the following
     * rule of thumb is that a lightweight post-processing shall be excuted by
     * default at each end of time steps and that an heavy post-processing (in
     * execution time and/or size of the generated file(s)) shall be executed at
     * post-processing times by default.
     */
    [[nodiscard]] virtual bool executePostProcessingTasks(
        Context &, const TimeStep &, const bool) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractPostProcessing() noexcept;
  };  // end of  class AbstractPostProcessing

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTPOSTPROCESSING_HXX */