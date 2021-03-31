/*!
 * \file   src/NewtonSolver.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#include <iomanip>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static void checkSolverOperator(
      const NonLinearEvolutionProblemImplementation<parallel> &p) {
    MFEM_ASSERT(p.Height() == p.Width(),
                "checkSolverOperator: "
                "a square operator is required.");
  }  // end of checkSolverOperator

#ifdef MFEM_USE_MPI

  template <>
  NewtonSolver<true>::NewtonSolver(
      NonLinearEvolutionProblemImplementation<true> &p)
      : IterativeSolver(p.getFiniteElementSpace().GetComm()), problem(p) {
    checkSolverOperator(p);
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions(
        [&p](const mfem::Vector &u) { return p.integrate(u); });
  }  // end of NewtonSolver

#endif /* MFEM_USE_MPI */

  template <>
  NewtonSolver<false>::NewtonSolver(
      NonLinearEvolutionProblemImplementation<false> &p)
      : problem(p) {
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions(
        [&p](const mfem::Vector &u) { return p.integrate(u); });
  }  // end of NewtonSolver

}  // end of namespace mfem_mgis
