/*!
 * \file   MFEMMGIS/AbstractConvergenceFailureHandler.hxx
 * \brief  This file declares the `AbstractConvergenceFailureHandler` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEMMGIS_ABSTRACTCONVERGENCEFAILUREHANDLER_HXX
#define LIB_MFEMMGIS_ABSTRACTCONVERGENCEFAILUREHANDLER_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief class used to determine a new time step in case
   * of convergence failure.
   */
  struct MFEM_MGIS_EXPORT AbstractConvergenceFailureHandler {
    /*!
     * \brief compute a new time increment after a convergence failure
     * \return the new time increment, empty on failure
     * \param[in, out] ctx: execution context
     * \param[in] dt: current time increment
     */
    virtual std::optional<real> getNewTimeIncrement(
        Context& ctx, const real dt) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractConvergenceFailureHandler();
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTCONVERGENCEFAILUREHANDLER_HXX */
