/*!
 * \file   src/SolverUtilities.cxx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#include "mfem/linalg/solvers.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/linalg/hypre.hpp"
#endif /* MFEM_USE_MPI */
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
  static void setSolverParametersImplementation(attributes::Throwing,
                                                SolverType& s,
                                                const Parameters& params) {
    using Problem = AbstractNonLinearEvolutionProblem;
    checkParameters(throwing, params, getIterativeSolverParametersList());
    if (contains(params, Problem::SolverVerbosityLevel)) {
      s.SetPrintLevel(
          get<int>(throwing, params, Problem::SolverVerbosityLevel));
    }
    if (contains(params, Problem::SolverRelativeTolerance)) {
      s.SetRelTol(
          get<double>(throwing, params, Problem::SolverRelativeTolerance));
    }
    if (contains(params, Problem::SolverAbsoluteTolerance)) {
      s.SetAbsTol(
          get<double>(throwing, params, Problem::SolverAbsoluteTolerance));
    }
    if (contains(params, Problem::SolverMaximumNumberOfIterations)) {
      s.SetMaxIter(
          get<int>(throwing, params, Problem::SolverMaximumNumberOfIterations));
    }
  }  // end of setSolverParametersImplementation

  bool setSolverParameters(Context& ctx,
                           IterativeSolver& s,
                           const Parameters& params) noexcept {
    try {
      setSolverParametersImplementation(throwing, s, params);
    } catch (...) {
      return mgis::registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of setSolverParameters

#ifdef MFEM_USE_PETSC
  bool setSolverParameters(Context& ctx,
                           mfem::PetscNonlinearSolver& s,
                           const Parameters& params) noexcept {
    try {
      setSolverParametersImplementation(throwing, s, params);
    } catch (...) {
      return mgis::registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }    // end of setSolverParameters
#endif /* MFEM_USE_PETSC */

  void setSolverParameters(attributes::Throwing,
                           IterativeSolver& s,
                           const Parameters& params) {
    setSolverParametersImplementation(throwing, s, params);
  }  // end of setSolverParameters

#ifdef MFEM_USE_PETSC
  void setSolverParameters(attributes::Throwing,
                           mfem::PetscNonlinearSolver& s,
                           const Parameters& params) {
    setSolverParametersImplementation(throwing, s, params);
  }    // end of setSolverParameters
#endif /* MFEM_USE_PETSC */

  bool hasConverged(const LinearSolver& ls) noexcept {
    if (const auto* isolver = dynamic_cast<const IterativeSolver*>(&ls);
        isolver != nullptr) {
      return isolver->GetConverged();
    }
#ifdef MFEM_USE_MPI
    auto check_hypre_solver_convergence = [](const auto* hptr) {
      int niter;
      hptr->GetNumIterations(niter);
      return niter < hptr->GetMaxIter();
    };
    if (const auto* hptr = dynamic_cast<const mfem::HyprePCG*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreGMRES*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
    if (const auto* hptr = dynamic_cast<const mfem::HypreFGMRES*>(&ls);
        hptr != nullptr) {
      return check_hypre_solver_convergence(hptr);
    }
#endif /* MFEM_USE_MPI */
    return true;
  }  // end of hasConverged

}  // end of namespace mfem_mgis
