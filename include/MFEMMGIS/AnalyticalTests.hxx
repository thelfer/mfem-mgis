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
   * \brief Compare the results to analytical solution with a specified threshold
   * \param[in] p: considered problem
   * \param[in] f: reference function to compare with
   * \param[in] params: set of parameters
   */
  MFEM_MGIS_EXPORT bool compareToAnalyticalSolution(
      NonLinearEvolutionProblem &p,
      std::function<void(mfem::Vector &, const mfem::Vector &)> f,
      const Parameters &params);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ANALYTICALTESTS_HXX */
