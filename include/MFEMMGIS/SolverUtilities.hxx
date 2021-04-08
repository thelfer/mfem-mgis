/*!
 * \file   MFEMMGIS/SolverUtilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_SOLVERUTILITIES_HXX
#define LIB_MFEM_MGIS_SOLVERUTILITIES_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief return the list of parameters for an interative solver
   * \note the following parameters are available:
   * - `AbstractNonLinearEvolutionProblem::SolverVerbosityLevel`, aka
   *   `"VerbosityLevel"`,
   * - `AbstractNonLinearEvolutionProblem::SolverRelativeTolerance`, aka
   *   `"RelativeTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverAbsoluteTolerance`, aka
   *   `"AbsoluteTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverMaximumNumberOfIterations`, aka
   *   `"MaximumNumberOfIterations"`,
   */
  MFEM_MGIS_EXPORT std::vector<std::string> getIterativeSolverParametersList();

  /*!
   * \brief set the parameters of an iterative solver
   * \param[in] s: iterative solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT void setSolverParameters(IterativeSolver&,
                                            const Parameters&);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SOLVERUTILITIES_HXX */
