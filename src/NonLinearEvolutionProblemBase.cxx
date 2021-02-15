/*!
 * \file   src/NonLinearEvolutionProblemBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemBase.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  NonLinearEvolutionProblemBase<true>::NonLinearEvolutionProblemBase(
      std::shared_ptr<FiniteElementDiscretization> fed)
      : NonLinearEvolutionProblemCommon(fed),
        mfem::ParNonlinearForm(&(fed->getFiniteElementSpace<true>())),
        solver(fed->getFiniteElementSpace<true>().GetComm()) {
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
  }  // end of NonLinearEvolutionProblemBase

  const FiniteElementSpace<true>&
  NonLinearEvolutionProblemBase<true>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<true>&
  NonLinearEvolutionProblemBase<true>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<true>();
  }  // end of getFiniteElementSpace

  NewtonSolver& NonLinearEvolutionProblemBase<true>::getSolver() {
    return this->solver;
  }  // end of getSolver

  void NonLinearEvolutionProblemBase<true>::solve(const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->setup();
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  NonLinearEvolutionProblemBase<true>::~NonLinearEvolutionProblemBase() = default;

#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblemBase<false>::NonLinearEvolutionProblemBase(
      std::shared_ptr<FiniteElementDiscretization> fed)
      : NonLinearEvolutionProblemCommon(fed),
        mfem::NonlinearForm(&(fed->getFiniteElementSpace<false>())) {
    this->solver.SetOperator(*(this));
    this->solver.iterative_mode = true;
  }  // end of NonLinearEvolutionProblemBase

  const FiniteElementSpace<false>&
  NonLinearEvolutionProblemBase<false>::getFiniteElementSpace() const {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  FiniteElementSpace<false>&
  NonLinearEvolutionProblemBase<false>::getFiniteElementSpace() {
    return this->fe_discretization->getFiniteElementSpace<false>();
  }  // end of getFiniteElementSpace

  NewtonSolver& NonLinearEvolutionProblemBase<false>::getSolver() {
    return this->solver;
  }  // end of getSolver

  void NonLinearEvolutionProblemBase<false>::solve(const real dt) {
    mfem::Vector zero;
    this->setTimeIncrement(dt);
    this->setup();
    this->solver.Mult(zero, this->u1);
    if (!this->solver.GetConverged()) {
      mgis::raise("Newton solver did not converge");
    }
  }  // end of solve

  NonLinearEvolutionProblemBase<false>::~NonLinearEvolutionProblemBase() = default;

}  // end of namespace mfem_mgis
