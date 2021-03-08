/*!
 * \file   include/MFEMMGIS/PostProcessing.hxx
 * \brief
 * \author Thomas Helfer
 * \date   08/03/2021
 */

#ifndef LIB_MFEM_MGIS_POSTPROCESSING_HXX
#define LIB_MFEM_MGIS_POSTPROCESSING_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblem;

  /*!
   * \brief base class for post-processing done at the end of each time step.
   */
  template <bool parallel>
  struct PostProcessing;

#ifdef MFEM_USE_MPI

  //! \brief partial specialisation for parallel post-processings
  template <>
  struct MFEM_MGIS_EXPORT PostProcessing<true> {
    /*!
     * \brief execute the post-processing
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void execute(NonLinearEvolutionProblem<true>&,
                         const real,
                         const real) = 0;
    //! \brief destructor
    virtual ~PostProcessing();
  };  // end of struct PostProcessing

#endif /* MFEM_USE_MPI */

  //! \brief partial specialisation for sequential post-processings
  template <>
  struct MFEM_MGIS_EXPORT PostProcessing<false> {
    /*!
     * \brief execute the post-processing
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void execute(NonLinearEvolutionProblem<false>&,
                         const real,
                         const real) = 0;
    //! \brief destructor
    virtual ~PostProcessing();
  };  // end of struct PostProcessing

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_POSTPROCESSING_HXX */
