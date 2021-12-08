/*!
 * \file   include/MFEMMGIS/AnalyticalTests.hxx
 * \brief
 * \author Thomas Helfer
 * \date   25/03/2021
 */

#ifndef LIB_MFEM_MGIS_ANALYTICALTESTS_HXX
#define LIB_MFEM_MGIS_ANALYTICALTESTS_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declaration
  struct NonLinearEvolutionProblem;

  /*!
   * \return the L2 norm of the difference between the unknowns at the end of
   * the time step and an analytical solution.
   * \param[in] p: considered problem
   * \param[in] f: reference function to compare with
   */
  MFEM_MGIS_EXPORT real computeL2ErrorAgainstAnalyticalSolution(
      NonLinearEvolutionProblem &,
      std::function<void(mfem::Vector &, const mfem::Vector &)>);

  /*!
   * \brief Compare the results to analytical solution with a specified
   * threshold using L2 norm.
   * \param[in] p: considered problem
   * \param[in] f: reference function to compare with
   * \param[in] params: set of parameters
   */
  MFEM_MGIS_EXPORT bool compareToAnalyticalSolution(
      NonLinearEvolutionProblem &,
      std::function<void(mfem::Vector &, const mfem::Vector &)>,
      const Parameters &);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ANALYTICALTESTS_HXX */
