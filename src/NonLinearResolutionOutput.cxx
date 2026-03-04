/*!
 * \file   MFEMMGIS/NonLinearResolutionOutput.cxx
 * \brief
 * \author Thomas Helfer
 * \date   04/03/2026
 */

#include "MFEMMGIS/NonLinearResolutionOutput.hxx"

namespace mfem_mgis {

  std::optional<Parameters> convertToParameters(
      const NonLinearResolutionOutput& output) noexcept {
    if (isInvalid(output)) {
      return {};
    }
    auto p = Parameters{};
    p.replaceOrInsert("initial_residual_norm", output.initial_residual_norm);
    p.replaceOrInsert("final_residual_norm", output.final_residual_norm);
    p.replaceOrInsert("iterations", output.iterations);
    return p;
  }  // end of convertToParameters

}  // end of namespace mfem_mgis
