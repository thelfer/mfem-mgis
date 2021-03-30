/*!
 * \file   src/NewtonSolver.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#include <iomanip>
#include <utility>
#include "MFEMMGIS/NewtonSolver.hxx"

namespace mfem_mgis {

  NewtonSolver::NewtonSolver() = default;

#ifdef MFEM_USE_MPI

  NewtonSolver::NewtonSolver(MPI_Comm c)
      : mfem::IterativeSolver(c) {}  // end of NewtonSolver

#endif /* MFEM_USE_MPI */

  void NewtonSolver::SetOperator(const mfem::Operator &op) {
    this->oper = &op;
    this->height = op.Height();
    this->width = op.Width();
    MFEM_ASSERT(height == width, "square Operator is required.");
  }  // end of SetOperator

  void NewtonSolver::setLinearSolver(LinearSolver &s) {
    this->prec = &s;
  }  // end of setLinearSolver

  void NewtonSolver::Mult(const mfem::Vector &b, mfem::Vector &x) const {
    MFEM_ASSERT(this->oper != nullptr, "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr, "the Solver is not set (use SetSolver).");

    mfem::Vector r, c;
    r.SetSize(this->oper->Width());
    c.SetSize(this->oper->Width());

    const auto have_b = (b.Size() == Height());
    this->final_iter = size_type{};
    this->final_norm = std::numeric_limits<real>::max();

    if (!this->processNewUnknownsEstimate(x)) {
      this->converged = 0;
      return;
    }

    if (!this->iterative_mode) {
      x = 0.0;
    }

    this->oper->Mult(x, r);
    if (have_b) {
      r -= b;
    }

    const auto norm0 = this->Norm(r);
    const auto norm_goal = std::max(rel_tol * norm0, abs_tol);

    const auto usesIterativeLinearSolver =
        dynamic_cast<const IterativeSolver *>(this->prec) != nullptr;
    this->prec->iterative_mode = false;

    // x_{i+1} = x_i - [DF(x_i)]^{-1} [F(x_i)-b]
    auto it = size_type{};
    auto norm = norm0;

    while (true) {
      MFEM_ASSERT(IsFinite(norm), "norm = " << norm);
      if (this->print_level >= 0) {
        mfem::out << "Newton iteration " << std::setw(2) << it
                  << " : ||r|| = " << norm;
        if (it > 0) {
          mfem::out << ", ||r||/||r_0|| = " << norm / norm0;
        }
        mfem::out << '\n';
      }
      Monitor(it, norm, r, x);

      if (norm <= norm_goal) {
        this->converged = 1;
        break;
      }

      if (it >= this->max_iter) {
        this->converged = 0;
        break;
      }

      this->prec->SetOperator(this->oper->GetGradient(x));
      this->prec->Mult(r, c);  // c = [DF(x_i)]^{-1} [F(x_i)-b]
      if (usesIterativeLinearSolver) {
        const auto &iprec =
            static_cast<const mfem::IterativeSolver &>(*(this->prec));
        if (!iprec.GetConverged()) {
          this->converged = 0;
          break;
        }
      }
      const double c_scale = ComputeScalingFactor(x, b);
      if (c_scale == 0.0) {
        this->converged = 0;
        break;
      }
      add(x, -c_scale, c, x);

      if (!this->processNewUnknownsEstimate(x)) {
        this->converged = 0;
        break;
      }

      this->oper->Mult(x, r);
      if (have_b) {
        r -= b;
      }
      norm = this->Norm(r);
      ++it;
    }
    this->final_iter = it;
    this->final_norm = norm;
  }  // end of Mult

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
