/*!
 * \file   src/Provider.cxx
 * \brief  This file implements the `Provider` class
 * \date   12/12/2022
 */

#include "MFEMMGIS/Provider.hxx"

namespace mfem_mgis {

  // ExitStatus Provider::reportInvalidResolveDependencyCall_(Context &ctx,
  // const ValueDependency &d)
  // {
  //   return ctx.registerErrorMessage("value dependency '" + d.getName() + "'
  //   on mesh set '" + d.getMeshSet().getName() +
  //                                   "' is not handled by this provider");
  // }    // end of reportInvalidResolveDependencyCall_
  //
  // ExitStatus Provider::reportInvalidResolveDependencyCall_(Context &ctx,
  // const NodalDependency &d)
  // {
  //   return ctx.registerErrorMessage("nodal dependency '" + d.getName() + "'
  //   on mesh set '" + d.getMeshSet().getName() +
  //                                   "' is not handled by this provider");
  // }    // end of reportInvalidResolveDependencyCall_
  //
  // ExitStatus Provider::reportInvalidResolveDependencyCall_(Context &ctx,
  // const IPDependency &d)
  // {
  //   return ctx.registerErrorMessage("dependency at integration point '" +
  //   d.getName() + "' on mesh set '" +
  //                                   d.getMeshSet().getName() + "' is not
  //                                   handled by this provider");
  // }    // end of reportInvalidResolveDependencyCall_

  Provider::~Provider() noexcept = default;

}  // end of namespace mfem_mgis
