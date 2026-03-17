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

  //! \brief the default time step computer
  struct MFEM_MGIS_EXPORT DefaultConvergenceFailureHandler
      : ConvergenceFailureHandlerBase {
    //! \brief constructor
    DefaultConvergenceFailureHandler() noexcept;
    //
    [[nodiscard]] std::optional<real> getNewTimeIncrement(
        Context &, const real) const noexcept override;
    //! \brief destructor
    ~DefaultConvergenceFailureHandler() override;
  };  // end of DefaultConvergenceFailureHandler

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_DEFAULTCONVERGENCEFAILUREHANDLER_HXX */
