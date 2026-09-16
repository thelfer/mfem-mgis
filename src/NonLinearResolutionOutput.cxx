/*!
 * \file   MFEMMGIS/NonLinearResolutionOutput.cxx
 * \brief
 * \author Thomas Helfer
 * \date   04/03/2026
 */

#include "MFEMMGIS/NonLinearResolutionOutput.hxx"

namespace mfem_mgis {

  std::optional<ComputeNextStateOutput> convertToComputeNextStateOutput(
      const NonLinearResolutionOutput& output) noexcept {
    if (isInvalid(output)) {
      return {};
    }
    auto solver = Parameters{};
    solver.replaceOrInsert("InitialResidualNorm", output.initial_residual_norm);
    solver.replaceOrInsert("FinalResidualNorm", output.final_residual_norm);
    solver.replaceOrInsert("NumberOfIterations", output.iterations);
    auto p = ComputeNextStateOutput{};
    p.replaceOrInsert("Solver", solver);
    return p;
  }  // end of convertToComputeNextStateOutput

}  // end of namespace mfem_mgis
