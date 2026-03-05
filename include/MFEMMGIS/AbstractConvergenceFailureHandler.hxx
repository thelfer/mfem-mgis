/*!
 * \file   MFEMMGIS/AbstractConvergenceFailureHandler.hxx
 * \brief  This class declares the `AbstractConvergenceFailureHandler` class
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
     * \return a new time increment on case of convergence failure
     * \param[in] ctx: execution context
     * \param[in] dt: current time incremnent
     */
    virtual std::optional<real> getNewTimeIncrement(
        Context &, const real) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractConvergenceFailureHandler();
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTCONVERGENCEFAILUREHANDLER_HXX */