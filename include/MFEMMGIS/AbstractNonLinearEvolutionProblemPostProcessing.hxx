/*!
 * \file   include/MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx
 * \brief  This file declares the
 * `AbstractNonLinearEvolutionProblemPostProcessing` class
 * \author Thomas Helfer
 * \date   08/03/2021
 */

#ifndef LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEMPOSTPROCESSING_HXX
#define LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEMPOSTPROCESSING_HXX

#include "MGIS/Context.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief base class for the post-processings of a non linear evolution
   * problem, executed at the initial time and at the end of each time step.
   */
  template <bool parallel>
  struct AbstractNonLinearEvolutionProblemPostProcessing;

#ifdef MFEM_USE_MPI

  //! \brief specialisation for parallel post-processings
  template <>
  struct MFEM_MGIS_EXPORT
      AbstractNonLinearEvolutionProblemPostProcessing<true> {
    /*!
     * \brief execute the post-processing at the initial time of the
     * simulation
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] virtual bool executeInitialPostProcessing(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<true>& p,
        const real t) noexcept = 0;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool execute(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<true>& p,
        const real t,
        const real dt) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearEvolutionProblemPostProcessing();
  };  // end of struct AbstractNonLinearEvolutionProblemPostProcessing

#endif /* MFEM_USE_MPI */

  //! \brief specialisation for sequential post-processings
  template <>
  struct MFEM_MGIS_EXPORT
      AbstractNonLinearEvolutionProblemPostProcessing<false> {
    /*!
     * \brief execute the post-processing at the initial time of the
     * simulation
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] virtual bool executeInitialPostProcessing(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<false>& p,
        const real t) noexcept = 0;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] virtual bool execute(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<false>& p,
        const real t,
        const real dt) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearEvolutionProblemPostProcessing();
  };  // end of struct AbstractNonLinearEvolutionProblemPostProcessing

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEMPOSTPROCESSING_HXX */
