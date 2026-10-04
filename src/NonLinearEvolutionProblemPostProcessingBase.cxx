/*!
 * \file src/NonLinearEvolutionProblemPostProcessingBase.cxx
 * \brief This file implements the `NonLinearEvolutionProblemPostProcessingBase`
 * class
 * \author Thomas Helfer
 * \date   04/10/2026
 */

#include "MFEMMGIS/PostProcessing/NonLinearEvolutionProblemPostProcessingBase.hxx"

namespace mfem_mgis {

  bool NonLinearEvolutionProblemPostProcessingBase<
      true>::hasExecuteInitialPostProcessingAlreadyBeenCalled() const noexcept {
    return this->executeInitialPostProcessingAlreadyCalled;
  }  // end of hasExecuteInitialPostProcessingAlreadyBeenCalled

  bool NonLinearEvolutionProblemPostProcessingBase<true>::
      executeInitialPostProcessing(
          Context&,
          NonLinearEvolutionProblemImplementation<true>&,
          const real) noexcept {
    this->executeInitialPostProcessingAlreadyCalled = true;
    return true;
  }  // end of executeInitialPostProcessing

  bool NonLinearEvolutionProblemPostProcessingBase<
      false>::hasExecuteInitialPostProcessingAlreadyBeenCalled() const noexcept {
    return this->executeInitialPostProcessingAlreadyCalled;
  }  // end of hasExecuteInitialPostProcessingAlreadyBeenCalled

  bool NonLinearEvolutionProblemPostProcessingBase<false>::
      executeInitialPostProcessing(
          Context&,
          NonLinearEvolutionProblemImplementation<false>&,
          const real) noexcept {
    this->executeInitialPostProcessingAlreadyCalled = true;
    return true;
  }  // end of executeInitialPostProcessing

}  // end of namespace mfem_mgis
