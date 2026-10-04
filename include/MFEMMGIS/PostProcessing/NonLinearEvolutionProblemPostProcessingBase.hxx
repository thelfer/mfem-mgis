/*!
 * \file MFEMMGIS/PostProcessing/NonLinearEvolutionProblemPostProcessingBase.hxx
 * \brief This file declares the `NonLinearEvolutionProblemPostProcessingBase`
 * class
 * \author Thomas Helfer
 * \date   04/10/2026
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_NONLINEAREVOLUTIONPROBLEMPOSTPROCESSINGBASE_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_NONLINEAREVOLUTIONPROBLEMPOSTPROCESSINGBASE_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief base class for most post-processings
   * of nonlinear evolution problems
   *
   * This class is mostly meant to ensure that `executeInitialPostProcessing` is
   * only called once.
   */
  template <bool parallel>
  struct NonLinearEvolutionProblemPostProcessingBase;

#ifdef MFEM_USE_MPI

  //! \brief partial specialization for parallel computations
  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemPostProcessingBase<true>
      : AbstractNonLinearEvolutionProblemPostProcessing<true> {
    bool hasExecuteInitialPostProcessingAlreadyBeenCalled()
        const noexcept override final;
    bool executeInitialPostProcessing(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<true>& p,
        const real t) noexcept override;

   private:
    /*!
     * \brief boolean stating if ` executeInitialPostProcessing` has been called
     */
    bool executeInitialPostProcessingAlreadyCalled = false;
  };  // end of NonLinearEvolutionProblemPostProcessingBase

#endif /* MFEM_USE_MPI */

  //! \brief partial specialization for sequential computations
  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemPostProcessingBase<false>
      : AbstractNonLinearEvolutionProblemPostProcessing<false> {
    bool hasExecuteInitialPostProcessingAlreadyBeenCalled()
        const noexcept override final;
    bool executeInitialPostProcessing(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<false>& p,
        const real t) noexcept override;

   private:
    /*!
     * \brief boolean stating if ` executeInitialPostProcessing` has been called
     */
    bool executeInitialPostProcessingAlreadyCalled = false;
  };  // end of NonLinearEvolutionProblemPostProcessingBase

}  // end of namespace mfem_mgis

#endif
