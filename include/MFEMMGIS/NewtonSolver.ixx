/*!
 * \file   MFEMMGIS/NewtonSolver.ixx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_NEWTONSOLVER_IXX
#define LIB_MFEM_MGIS_NEWTONSOLVER_IXX

#include <utility>
#include "MGIS/Raise.hxx"

namespace mfem_mgis {

  template <bool parallel>
  void NewtonSolver<parallel>::SetOperator(const mfem::Operator &) {
    mgis::raise("NewtonSolver<parallel>::SetOperator: invalid call");
  }  // end of SetOperator

  template <bool parallel>
  void NewtonSolver<parallel>::SetPreconditioner(Solver &) {
    mgis::raise("NewtonSolver<parallel>::SetOperator: invalid call");
  }  // end of SetPreconditioner

  template <bool parallel>
  void NewtonSolver<parallel>::Mult(const mfem::Vector &,
                                    mfem::Vector &) const {
    mgis::raise("NewtonSolver<parallel>::Mult: invalid call");
  }  // end of Mult

  template <bool parallel>
  void NewtonSolver<parallel>::setLinearSolver(
      std::unique_ptr<LinearSolver> s) {
    this->linear_solver = std::move(s);
    this->prec = this->linear_solver.get();
  }  // end of setLinearSolver

  template <bool parallel>
  void NewtonSolver<parallel>::solve() const {
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");

    auto &x = this->problem.getUnknownsAtEndOfTheTimeStep();

    mfem::Vector r, c;
    r.SetSize(this->oper->Width());
    c.SetSize(this->oper->Width());

    const auto usesIterativeLinearSolver =
        dynamic_cast<const IterativeSolver *>(this->prec) != nullptr;
    this->prec->iterative_mode = false;

    auto updateResidual = [this, &r, &x] {
      this->oper->Mult(x, r);
      return this->Norm(r);
    };

    auto computeNewtonCorrection = [this, &c, &r, &x,
                                    usesIterativeLinearSolver] {
      this->prec->SetOperator(this->oper->GetGradient(x));
      this->prec->Mult(r, c);  // c = [DF(x_i)]^{-1} [F(x_i)-b]
      if (usesIterativeLinearSolver) {
        const auto &iprec =
            static_cast<const mfem::IterativeSolver &>(*(this->prec));
        if (!iprec.GetConverged()) {
          return false;
        }
      }
      return true;
    };

    this->final_iter = size_type{};
    this->final_norm = std::numeric_limits<real>::max();

    if (!this->processNewUnknownsEstimate(x)) {
      this->converged = 0;
      return;
    }
    const auto norm0 = updateResidual();

    if (this->prediction) {
    }
    //
    //       if (!computeNewtonCorrection()) {
    //         this->converged = false;
    //         return;
    //       }
    //       add(x, -1, c, x);
    //     } else {
    //       if (!this->processNewUnknownsEstimate(x)) {
    //         this->converged = 0;
    //         return;
    //       }
    //       norm0 = updateResidual();
    //     }

    //
    const auto norm_goal = std::max(rel_tol * norm0, abs_tol);
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
      if (!computeNewtonCorrection()) {
        this->converged = 0;
        break;
      }
      //
      // x_{i+1} = x_i - c * [DF(x_i)]^{-1} [F(x_i)-b]
      add(x, -1, c, x);

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

  template <bool parallel>
  void NewtonSolver<parallel>::addNewUnknownsEstimateActions(
      std::function<bool(const mfem::Vector &)> a) {
    this->nue_actions.push_back(std::move(a));
  }  // end of addNewUnknownsEstimateActions

  template <bool parallel>
  bool NewtonSolver<parallel>::processNewUnknownsEstimate(const mfem::Vector &u) const {
    for (const auto &a : this->nue_actions) {
      if (!a(u)) {
        return false;
      }
    }
    return true;
  }  // end of processNewUnknownsEstimate

  template <bool parallel>
  NewtonSolver<parallel>::~NewtonSolver() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NEWTONSOLVER_IXX */
