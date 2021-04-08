/*!
 * \file   include/MFEMMGIS/ComputeStressStrain.hxx
 * \brief
 * \author Guillaume Latu
 * \date   07/04/2021
 */

#ifndef LIB_COMPUTESTRESSSTRAIN_HXX
#define LIB_COMPUTESTRESSSTRAIN_HXX

#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing computing the resultant of the inner forces on a
   * boundary.
   *
   * The following parameters are required:
   *
   * - `Boundary`: identifier of the boundary
   * - `OutputFileName`: name of the output file
   */
  template <bool parallel>
  struct ComputeStressStrain;

  /*!
   * \brief a base class for the `ComputeStressStrain`
   * post-processing
   */
  struct ComputeStressStrainCommon {
    /*!
     * \brief constructor
     * \param[in] i: component to be exported
     */
    ComputeStressStrainCommon(const size_type, std::shared_ptr<FiniteElementDiscretization>);

    //! \brief component to be exported
    const size_type icomp;
    //! \brief underlying finite element discretization
    const std::shared_ptr<FiniteElementDiscretization> fed;
  };  // end of struct ComputeStressStrainCommon
  
#ifdef MFEM_USE_MPI

  /*!
   * \brief partial specialisation of the `ComputeStressStrain`
   * post-processing in parallel
   */
  template <>
  struct MFEM_MGIS_EXPORT ComputeStressStrain<true> final
    : public PostProcessing<true>,
      protected ComputeStressStrainCommon {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ComputeStressStrain(
        NonLinearEvolutionProblemImplementation<true>&, const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<true>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ComputeStressStrain() override;
  };  // end of struct ComputeStressStrain

#endif /* MFEM_USE_MPI */

  /*!
   * \brief partial specialisation of the `ComputeStressStrain`
   * post-processing in sequential
   */
  template <>
  struct MFEM_MGIS_EXPORT ComputeStressStrain<false> final
      : public PostProcessing<false>,
        protected ComputeStressStrainCommon {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ComputeStressStrain(
        NonLinearEvolutionProblemImplementation<false>&, const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<false>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ComputeStressStrain() override;
  };  // end of struct ComputeStressStrain

}  // end of namespace mfem_mgis

#endif /* LIB_COMPUTESTRESSSTRAIN_HXX */
