/*!
 * \file   MFEMMGIS/AbstractSimulationMonitor.hxx
 * \brief  This file declares the AbstractSimulationMonitor class
 * \date   11/09/2024
 */

#ifndef LIB_MFEMMGIS_ABSTRACT_SIMULATION_MONITOR_HXX
#define LIB_MFEMMGIS_ABSTRACT_SIMULATION_MONITOR_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declarations
  struct SimulationOutput;

  /*!
   * \brief an abstract class used to monitor a simulation at the end of each
   * successful time step or at the end of the simulation
   *
   * Such monitors can be used to report meaningful information to the user
   * about the convergence or the resources usage.
   */
  struct MFEM_MGIS_EXPORT AbstractSimulationMonitor {
    /*!
     * \brief execute the monitor
     *
     * \param[in, out] ctx: execution context
     * \param[in] output: outputs of the simulation
     * \param[in] endOfSimulation: boolean
     * stating if the simulation is over
     * \return true on success
     */
    virtual bool execute(Context& ctx,
                         const SimulationOutput& output,
                         const bool endOfSimulation) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractSimulationMonitor() noexcept;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACT_SIMULATION_MONITOR_HXX */
