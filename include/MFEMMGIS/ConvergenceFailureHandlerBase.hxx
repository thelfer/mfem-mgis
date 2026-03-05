/*!
 * \file   MFEMMGIS/ConvergenceFailureHandlerBase.hxx
 * \brief  This class declares the `ConvergenceFailureHandlerBase` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEM_MGIS_CONVERGENCEFAILUREHANDLERBASE_HXX
#define LIB_MFEM_MGIS_CONVERGENCEFAILUREHANDLERBASE_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractConvergenceFailureHandler.hxx"

namespace mfem_mgis {

  //! \brief a common class for most time step computers
  struct MFEM_MGIS_EXPORT ConvergenceFailureHandlerBase
      : AbstractConvergenceFailureHandler {
    //!\brief constructor
    ConvergenceFailureHandlerBase() noexcept;
    //! \brief destructor
    ~ConvergenceFailureHandlerBase() noexcept override;
  };

}  // end of namespace mfem_mgis

#endif /*LIB_MFEM_MGIS_CONVERGENCEFAILUREHANDLERBASE_HXX */