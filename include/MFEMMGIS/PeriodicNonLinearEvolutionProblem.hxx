/*!
 * \file   include/MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct PeriodicNonLinearEvolutionProblem;

#ifdef MFEM_USE_MPI

  template <>
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem<true>
      : public NonLinearEvolutionProblem<true> {
    /*!
     * \brief an helper function fix the degree of freedom of a point
     * \param[in] p: non linear evolution problem
     */
    static void setBoundaryConditions(
        mfem_mgis::NonLinearEvolutionProblemBase<true>&);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>);
    //! \brief destructor
    ~PeriodicNonLinearEvolutionProblem() override;
  };  // end of struct PeriodicNonLinearEvolutionProblem

#endif /* MFEM_USE_MPI */

  template <>
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem<false>
      : public NonLinearEvolutionProblem<false> {
    /*!
     * \brief an helper function fix the degree of freedom of a point
     * \param[in] p: non linear evolution problem
     */
    static void setBoundaryConditions(
        mfem_mgis::NonLinearEvolutionProblemBase<false>&);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>);
    //! \brief destructor
    ~PeriodicNonLinearEvolutionProblem() override;
  };  // end of struct PeriodicNonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM */
