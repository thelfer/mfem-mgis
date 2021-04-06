/*!
 * \file   src/SolverUtilities.cxx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#include "mfem/linalg/solvers.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"

namespace mfem_mgis {

  void setSolverParameters(IterativeSolver& solver, const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    checkParameters(params, {Problem::SolverVerbosityLevel,
                             Problem::SolverRelativeTolerance,
                             Problem::SolverAbsoluteTolerance,
                             Problem::SolverMaximumNumberOfIterations});
    if (contains(params, Problem::SolverVerbosityLevel)) {
      solver.SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
    }
    if (contains(params, Problem::SolverRelativeTolerance)) {
      solver.SetRelTol(get<double>(params, Problem::SolverRelativeTolerance));
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      solver.SetAbsTol(get<double>(params, Problem::SolverAbsoluteTolerance));
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      solver.SetMaxIter(
          get<int>(params, Problem::SolverMaximumNumberOfIterations));
    }
  }  // end of setSolverParameters

}  // end of namespace mfem_mgis