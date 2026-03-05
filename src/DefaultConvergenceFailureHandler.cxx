/*!
 * \file   MFEMMGIS/DefaultConvergenceFailureHandler.cxx
 * \brief  This file implements the `DefaultConvergenceFailureHandler` class
 * \date   08/12/2023
 */

#include "MFEMMGIS/DefaultConvergenceFailureHandler.hxx"

namespace mfem_mgis {

  DefaultConvergenceFailureHandler::
      DefaultConvergenceFailureHandler() noexcept = default;

  std::optional<real> DefaultConvergenceFailureHandler::getNewTimeIncrement(
      Context &, const real dt) const noexcept {
    return dt / 2;
  }  // end of getNextTimeStep

  DefaultConvergenceFailureHandler::~DefaultConvergenceFailureHandler() =
      default;

}  // end of namespace mfem_mgis
