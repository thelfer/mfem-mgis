/*!
 * \file   MFEMMGIS/AbstractTimeIncrementComputer.hxx
 * \brief  This class declares the `AbstractTimeIncrementComputer` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEMMGIS_ABSTRACTTIMEINCREMENTCOMPUTER_HXX
#define LIB_MFEMMGIS_ABSTRACTTIMEINCREMENTCOMPUTER_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief class used to determine a new time increment at the
   * beginning of a new time step.
   */
  struct MFEM_MGIS_EXPORT AbstractTimeIncrementComputer {
    /*!
     * \brief method called before the start of a simulation
     * \param[in] ctx: execution context
     */
    [[nodiscard]] virtual bool initialize(Context &) noexcept = 0;
    /*!
     * \brief this method is called at the end of a time step,
     * before updating the state of the system.
     *
     * This method allows the time step computer to predict
     * the next time increment.
     *
     * \param[in] ctx: execution context
     */
    [[nodiscard]] virtual bool prepareNextTimeStep(Context &) noexcept = 0;
    /*!
     * \return the next time increment
     * \param[in] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    [[nodiscard]] virtual std::optional<real> getNextTimeIncrement(
        Context &, const real, const real) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractTimeIncrementComputer();
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTTIMEINCREMENTCOMPUTER_HXX */
