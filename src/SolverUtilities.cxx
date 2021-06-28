/*!
 * \file   src/SolverUtilities.cxx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#include "mfem/linalg/solvers.hpp"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/SolverUtilities.hxx"

namespace mfem_mgis {

  std::vector<std::string> getIterativeSolverParametersList() {
    using Problem = AbstractNonLinearEvolutionProblem;
    return {Problem::SolverVerbosityLevel, Problem::SolverRelativeTolerance,
            Problem::SolverAbsoluteTolerance,
            Problem::SolverMaximumNumberOfIterations};
  }  // end of getIterativeSolverParametersList

  template <typename SolverType>
  static void setSolverParametersImplementation(SolverType& s,
                                                const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    checkParameters(params, getIterativeSolverParametersList());
    if (contains(params, Problem::SolverVerbosityLevel)) {
      s.SetPrintLevel(get<int>(params, Problem::SolverVerbosityLevel));
    }
    if (contains(params, Problem::SolverRelativeTolerance)) {
      s.SetRelTol(get<double>(params, Problem::SolverRelativeTolerance));
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      s.SetAbsTol(get<double>(params, Problem::SolverAbsoluteTolerance));
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      s.SetMaxIter(get<int>(params, Problem::SolverMaximumNumberOfIterations));
    }
  }  // end of setSolverParametersImplementation

  void setSolverParameters(IterativeSolver& s, const Parameters& params) {
    setSolverParametersImplementation(s, params);
  }  // end of setSolverParameters

#ifdef MFEM_USE_PETSC
  void setSolverParameters(mfem::PetscNonlinearSolver& s,
                           const Parameters& params) {
    setSolverParametersImplementation(s, params);
  }    // end of setSolverParameters
#endif /* MFEM_USE_PETSC */

}  // end of namespace mfem_mgis