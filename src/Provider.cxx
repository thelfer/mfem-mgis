/*!
 * \file   src/Provider.cxx
 * \brief  This file implements the `Provider` class
 * \author Thomas Helfer
 * \date   12/12/2022
 */

#include <string>
#include "MFEMMGIS/Provider.hxx"
#include "MFEMMGIS/Dependency.hxx"

namespace mfem_mgis {

  bool Provider::reportInvalidResolveDependencyCall(
      Context &ctx, const QPDependency &d) noexcept {
    return ctx.registerErrorMessage("dependency at integration point '" +
                                    d.getName() + "' on material '" +
                                    std::to_string(d.getMaterialIdentifier()) +
                                    "' is not handled by this provider");
  }  // end of reportInvalidResolveDependencyCall

  Provider::~Provider() noexcept = default;

}  // end of namespace mfem_mgis
