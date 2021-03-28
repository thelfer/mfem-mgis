/*!
 * \file   src/Config.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/02/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  [[noreturn]] void reportUnsupportedParallelComputations() {
    mgis::raise(
        "FiniteElementDiscretization::reportUnsupportedParallelComputations: "
        "unsupported parallel computations");
  }

}  // namespace mfem_mgis
