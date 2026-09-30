/*!
 * \file   MFEMMGIS/DefaultConvergenceFailureHandler.hxx
 * \brief  This file declares the `DefaultConvergenceFailureHandler` class
 * \date   08/12/2023
 */

#ifndef LIB_MFEM_MGIS_DEFAULTCONVERGENCEFAILUREHANDLER_HXX
#define LIB_MFEM_MGIS_DEFAULTCONVERGENCEFAILUREHANDLER_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ConvergenceFailureHandlerBase.hxx"

namespace mfem_mgis {

  //! \brief the default convergence failure handler
  struct MFEM_MGIS_EXPORT DefaultConvergenceFailureHandler
      : ConvergenceFailureHandlerBase {
    //! \brief constructor
    DefaultConvergenceFailureHandler() noexcept;
    /*!
     * \brief divide the current time increment by two
     * \return half of the current time increment
     * \param[in, out] ctx: execution context
     * \param[in] dt: current time increment
     */
    [[nodiscard]] std::optional<real> getNewTimeIncrement(
        Context& ctx, const real dt) const noexcept override;
    //! \brief destructor
    ~DefaultConvergenceFailureHandler() override;
  };  // end of DefaultConvergenceFailureHandler

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_DEFAULTCONVERGENCEFAILUREHANDLER_HXX */
