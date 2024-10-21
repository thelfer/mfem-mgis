/*!
 * \file   src/NewtonSolver.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#include <iomanip>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/NewtonSolver.hxx"
#include "MFEMMGIS/Profiler.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static void checkSolverOperator(
      const NonLinearEvolutionProblemImplementation<parallel> &p) {
    MFEM_ASSERT(p.Height() == p.Width(),
                "checkSolverOperator: "
                "a square operator is required.");
    static_cast<void>(p);
  }  // end of checkSolverOperator

#ifdef MFEM_USE_MPI

  NewtonSolver::NewtonSolver(NonLinearEvolutionProblemImplementation<true> &p)
      : IterativeSolver(p.getFiniteElementSpace().GetComm()) {
    checkSolverOperator(p);
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions([&p](const mfem::Vector &u) {
      return p.integrate(
          u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);
    });
  }  // end of NewtonSolver

#endif /* MFEM_USE_MPI */

  NewtonSolver::NewtonSolver(
      NonLinearEvolutionProblemImplementation<false> &p) {
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions([&p](const mfem::Vector &u) {
      return p.integrate(
          u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);
    });
  }  // end of NewtonSolver

  void NewtonSolver::SetOperator(const mfem::Operator &) {
    raise("NewtonSolver::SetOperator: invalid call");
  }  // end of SetOperator

  void NewtonSolver::SetPreconditioner(Solver &) {
    raise("NewtonSolver::SetOperator: invalid call");
  }  // end of SetPreconditioner

  void NewtonSolver::setLinearSolver(LinearSolver &s) {
    this->prec = &s;
    this->prec->iterative_mode = false;
  }  // end of setLinearSolver

  real NewtonSolver::GetInitialNorm() const {
    return this->initial_norm;
  }  // end of GetInitialNorm

  void NewtonSolver::Mult(const mfem::Vector &, mfem::Vector &x) const {
    CatchTimeSection("NewtonSolver::Mult");
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");

    mfem::Vector r, c;
    r.SetSize(this->oper->Width());
    c.SetSize(this->oper->Width());

    auto updateResidual = [this, &r, &x] {
      this->computeResidual(r, x);
      return this->Norm(r);
    };

    this->final_iter = size_type{};
    this->final_norm = std::numeric_limits<real>::max();

    if (!this->processNewUnknownsEstimate(x)) {
      this->converged = 0;
      return;
    }
    this->initial_norm = updateResidual();

    //
    const auto norm_goal = std::max(rel_tol * (this->initial_norm), abs_tol);
    auto it = size_type{};
    auto norm = this->initial_norm;

    while (true) {
      MFEM_ASSERT(mfem::IsFinite(norm), "norm = " << norm);
      if (this->print_level >= 0) {
        if (mfem_mgis::getMPIrank() == 0)
          mfem::out << "Newton iteration " << std::setw(2) << it
                    << " : ||r|| = " << norm;
        if (it > 0) {
          if (mfem_mgis::getMPIrank() == 0)
            mfem::out << ", ||r||/||r_0|| = " << norm / (this->initial_norm);
        }
        if (mfem_mgis::getMPIrank() == 0) mfem::out << '\n';
      }
      this->Monitor(it, norm, r, x);
      //
      if (norm <= norm_goal) {
        this->converged = 1;
        break;
      }
      //
      if (it >= this->max_iter) {
        this->converged = 0;
        break;
      }
      //
      if (!this->computeNewtonCorrection(c, r, x)) {
        this->converged = 0;
        break;
      }
      //
      // x_{i+1} = x_i - c * [DF(x_i)]^{-1} [F(x_i)-b]
      //      add(x, -1, c, x);
      x -= c;

      if (!this->processNewUnknownsEstimate(x)) {
        this->converged = 0;
        break;
      }

      updateResidual();
      norm = this->Norm(r);
      ++it;
    }
    this->final_iter = it;
    this->final_norm = norm;
  }  // end of Mult

  void NewtonSolver::computeResidual(mfem::Vector &r,
                                     const mfem::Vector &u) const {
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    this->oper->Mult(u, r);
  }  // end of NewtonSolver::computeResidual

  bool NewtonSolver::computeNewtonCorrection(mfem::Vector &c,
                                             const mfem::Vector &r,
                                             const mfem::Vector &u) const {
    CatchTimeSection("NewtonSolver::computeNewtonCorrection");
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    const auto usesIterativeLinearSolver =
        dynamic_cast<const IterativeSolver *>(this->prec) != nullptr;
    this->prec->SetOperator(this->getJacobian(u));
    this->prec->Mult(r, c);  // c = [DF(x_i)]^{-1} [F(x_i)-b]
    if (usesIterativeLinearSolver) {
      const auto &iprec =
          static_cast<const mfem::IterativeSolver &>(*(this->prec));
      return iprec.GetConverged();
    }
    return true;
  }  // end of computeNewtonCorrection

  mfem::Operator &NewtonSolver::getJacobian(const mfem::Vector &u) const {
    CatchTimeSection("NewtonSolver::getJacobian");
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    return this->oper->GetGradient(u);
  }  // end of getJacobian

  void NewtonSolver::addNewUnknownsEstimateActions(
      std::function<bool(const mfem::Vector &)> a) {
    this->nue_actions.push_back(std::move(a));
  }  // end of addNewUnknownsEstimateActions

  bool NewtonSolver::processNewUnknownsEstimate(const mfem::Vector &u) const {
    for (const auto &a : this->nue_actions) {
      if (!a(u)) {
        return false;
      }
    }
    return true;
  }  // end of processNewUnknownsEstimate

  NewtonSolver::~NewtonSolver() = default;

}  // end of namespace mfem_mgis
