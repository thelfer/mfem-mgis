/*!
 * \file   src/AbstractNonLinearSolver.cxx
 * \brief
 * \author Thomas Helfer
 * \date   20/09/2026
 */

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearSolvers/AbstractNonLinearSolver.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  AbstractNonLinearSolver::AbstractNonLinearSolver(
      NonLinearEvolutionProblemImplementation<true> &p)
      : IterativeSolver(p.getFiniteElementSpace().GetComm()) {
  }  // end of AbstractNonLinearSolver

#endif /* MFEM_USE_MPI */

  AbstractNonLinearSolver::AbstractNonLinearSolver(
      NonLinearEvolutionProblemImplementation<false> &) {
  }  // end of AbstractNonLinearSolver

  AbstractNonLinearSolver::~AbstractNonLinearSolver() = default;

}  // end of namespace mfem_mgis
