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

  MFEM_MGIS_EXPORT bool compareToAnalyticalSolution(
      NonLinearEvolutionProblem &,
      std::function<void(mfem::Vector &, const mfem::Vector &)>,
      const Parameters &);

  //   MFEM_MGIS_EXPORT bool compareToAnalyticalSolution(
  //       NonLinearEvolutionProblem &,
  //       std::function<void(mfem::Vector &, const mfem::Vector &, const
  //       double)>, const Parameters);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ANALYTICALTESTS_HXX */
