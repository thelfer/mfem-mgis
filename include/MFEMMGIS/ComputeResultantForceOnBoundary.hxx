/*!
 * \file   include/MFEMMGIS/ComputeResultantForceOnBoundary.hxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_COMPUTERESULTANTFORCEONBOUNDARY_HXX
#define LIB_COMPUTERESULTANTFORCEONBOUNDARY_HXX

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
  struct ComputeResultantForceOnBoundary;

  /*!
   * \brief a base class for the `ComputeResultantForceOnBoundary`
   * post-processing
   */
  struct ComputeResultantForceOnBoundaryCommon {
    /*!
     * \brief constructor
     * \param[in] edofs: elements degrees of freedom
     * \param[in] i: boundary identifier
     */
    ComputeResultantForceOnBoundaryCommon(
        std::vector<std::pair<size_type, std::vector<std::vector<size_type>>>>,
        const size_type);
    /*!
     * \brief return a structure which associates the global number of the
     * selected elements to the local indexes of its degrees of freedom sorted
     * by components.
     */
    const std::vector<std::pair<size_type,  // element number
                                std::vector<std::vector<size_type>>>>
        elts_dofs;
    //! \brief boundary identifier
    const size_type bid;
    //! \brief output file
    std::ofstream out;
  };  // end of struct ComputeResultantForceOnBoundaryCommon

#ifdef MFEM_USE_MPI

  /*!
   * \brief partial specialisation of the `ComputeResultantForceOnBoundary`
   * post-processing in parallel
   */
  template <>
  struct MFEM_MGIS_EXPORT ComputeResultantForceOnBoundary<true> final
      : public PostProcessing<true>,
        protected ComputeResultantForceOnBoundaryCommon {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ComputeResultantForceOnBoundary(
        NonLinearEvolutionProblemImplementation<true>&, const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<true>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ComputeResultantForceOnBoundary() override;
  };  // end of struct ComputeResultantForceOnBoundary

#endif /* MFEM_USE_MPI */

  /*!
   * \brief partial specialisation of the `ComputeResultantForceOnBoundary`
   * post-processing in sequential
   */
  template <>
  struct MFEM_MGIS_EXPORT ComputeResultantForceOnBoundary<false> final
      : public PostProcessing<false>,
        protected ComputeResultantForceOnBoundaryCommon {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ComputeResultantForceOnBoundary(
        NonLinearEvolutionProblemImplementation<false>&, const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<false>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ComputeResultantForceOnBoundary() override;
  };  // end of struct ComputeResultantForceOnBoundary

}  // end of namespace mfem_mgis

#endif /* LIB_COMPUTERESULTANTFORCEONBOUNDARY_HXX */
