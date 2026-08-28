/*!
 * \file   src/NewtonSolver.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#include <array>
#include <iomanip>
#include <utility>
#include "MGIS/Raise.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/IntegrationType.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"
#include "MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx"

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
      : NonLinearSolverBase(p) {}  // end of NewtonSolver

#endif /* MFEM_USE_MPI */

  NewtonSolver::NewtonSolver(NonLinearEvolutionProblemImplementation<false> &p)
      : NonLinearSolverBase(p) {}  // end of NewtonSolver

  
  void NewtonSolver::addAdditionalConvergenceCheck(std::shared_ptr<nonlinear_solver::AbstractAdditionalConvergenceCriterion> cv_check) {
    // TODO
    //CatchTimeSection("NS::addAdditionalConvergenceCheck);
    this->acc_actions.push_back(std::move(cv_check));

  } // end of addAdditionalConvergenceCheck

  /* // Unused ? Because of the change from std::function to a struct */
  std::optional<bool> NewtonSolver::processAdditionalConvergenceCheck(Context& ctx, const nonlinear_solver::AbstractAdditionalConvergenceCriterion::CheckArguments& s) const  {
    //TODO
    //CatchTimeSection("NS::processAdditionalConvergenceCheck");
    bool cv = s.converged;
    // a->check must be called (for each element of the list, in case it manipulates some values as a side effect)
      for (auto& a : this->acc_actions) {
          std::optional<bool> result = a->check(ctx,s);
          if (isInvalid(result)){
              return {};
          }
          cv = cv && *result; 
      }
    return cv;
  }  // end of processAdditionalConvergenceCheck

  void NewtonSolver::processAdditionalConvergenceReset()  {
    //TODO
    //CatchTimeSection("NS::processAdditionalConvergenceReset");
    for (auto& a : this->acc_actions) {
      a->reset();
    }
  }  // end of processAdditionalConvergenceReset
  
  void NewtonSolver::processAdditionalConvergenceHelper()  {
    //TODO
    //CatchTimeSection("NS::processAdditionalConvergenceHelper");
    for (auto& a : this->acc_actions) {
      a->helper();
    }
  }  // end of processAdditionalConvergenceReset

  void NewtonSolver::Mult(const mfem::Vector &, mfem::Vector &x) const {
    auto profiler_mult =
        this->ctx_ptr != nullptr
            ? this->ctx_ptr->startNewProfiling(
                  "NS::Mult", this->ctx_ptr->isProfilingEnabled())
            : mgis::ProfilingSection{};  // fallback RAII
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    // log stream
    auto &log = [this]() -> std::ostream & {
      if (this->ctx_ptr == nullptr) {
        return getDefaultLogStream();
      }
      return this->ctx_ptr->log();
    }();
    // boolean stating if messages shall be displayed
    auto shall_print = [this]() -> bool {
      if (this->print_level >= 0) {
#ifdef MFEM_USE_MPI
        // We must check that a communicator has been set.
        // This is not the case in sequential computations
        if (this->GetComm() != MPI_COMM_NULL) {
          int rank = 0;
          MPI_Comm_rank(this->GetComm(), &rank);
          return rank == 0;
        }
        return true;
#else  /* MFEM_USE_MPI */
        return true;
#endif /* MFEM_USE_MPI */
      }
      if (this->ctx_ptr != nullptr) {
        return this->ctx_ptr->getVerbosityLevel() >=
               VerbosityLevel::verboseLevel3;
      }
      return false;
    }();

    mfem::Vector r;  // residual vector
    mfem::Vector c;  // opposite of the Newton's correction
    r.SetSize(this->oper->Width());
    c.SetSize(this->oper->Width());

    auto updateResidual = [this, &r, &x] {
      if (this->ctx_ptr != nullptr) {
        CatchTimeSection(*(this->ctx_ptr), "NS::computeResidual");
        this->computeResidual(r, x);
      } else {
        this->computeResidual(r, x);
      }
      return this->Norm(r);
    };

    this->final_iter = size_type{};
    this->final_norm = std::numeric_limits<real>::max();

    if (!this->processNewUnknownsEstimate(x)) {
      this->converged = 0;
      return;
    }
    auto norm = updateResidual();
    if (!this->reference_residual_norm.has_value()) {
      this->reference_residual_norm = norm;
    }
    // this data member is not used, but we define it by
    // consistency
    this->initial_norm = *(this->reference_residual_norm);
    auto previous_norms =
        std::array<real, 2u>{this->initial_norm, this->initial_norm};

    const auto norm_goal =
        std::max(rel_tol * (*(this->reference_residual_norm)), abs_tol);
    auto it = size_type{};

    this->converged = 0;
    while (true) {
      auto profiler_while =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "NS::Mult::WhileLoop", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      MFEM_ASSERT(mfem::IsFinite(norm), "norm = " << norm);
      if (shall_print) {
        log << "Newton iteration " << std::setw(2) << it
            << " : ||r|| = " << norm
            << ", ||r||/||r_0|| = " << norm / (*(this->reference_residual_norm))
            << '\n';
      }
      this->Monitor(it, norm, r, x);
      //
      auto result = this->processAdditionalConvergenceCheck(*this->ctx_ptr, {
          .residual_norm = norm,
          .reference_residual_norm = this->reference_residual_norm.value(),
          .iter = it ,
          .max_iter = this->max_iter,
          .converged = norm <= norm_goal,
          .u = x
          }
          );     
      if (isInvalid(result)){
        this->converged=false;
        break;
      }
      this->converged = *result;
      if (this->converged){
          break;
      }
      //
      if (it >= this->max_iter) {
        break;
      }
      //
      if (!this->computeNewtonCorrection(c, r, x)) {
        break;
      }
      //
      // x_{i+1} = x_i - c * [DF(x_i)]^{-1} [F(x_i)-b]
      //      add(x, -1, c, x);
      x -= c;

      if (!this->processNewUnknownsEstimate(x)) {
        ++it;
        // basic line-search in case of integration failure
        //
        // we keep the current direction, but reduce the amplitude by a factor
        // two until we find an estimate of the solution that do not lead
        // to an integration failure or that the number of iterations reaches
        // the maximum value
        while (true) {
          if (it >= this->max_iter) {
            break;
          }
          if (shall_print) {
            log << "Newton iteration " << std::setw(2) << it
                << ": reducing the amplitude of the correction by a "
                   "factor 2\n";
          }
          c *= real{1} / 2;
          x += c;
          if (this->processNewUnknownsEstimate(x)) {
            break;
          }
          ++it;
        }
        if (it >= this->max_iter) {
          break;
        }
      }
     
      updateResidual();
      previous_norms[0] = previous_norms[1];
      previous_norms[1] = norm;
      norm = this->Norm(r);

     ++it;
    }
    this->final_iter = it;
    this->final_norm = norm;
    if (this->converged == 1) {
      // estimation of the convergence order
      if (shall_print) {
        if ((it >= 2) && (norm > 100 * std::numeric_limits<real>::min()) &&
            (previous_norms[0] > 100 * std::numeric_limits<real>::min()) &&
            (previous_norms[1] > 100 * std::numeric_limits<real>::min())) {
          const auto e1 = std::log(norm / previous_norms[1]);
          const auto e2 = std::log(previous_norms[1] / previous_norms[0]);
          if (std::abs(e2) > 100 * std::numeric_limits<real>::min()) {
            log << "Convergence order " << e1 / e2 << "\n\n";
          } else {
            log << "Convergence order undefined\n\n";
          }
        } else {
          log << "Convergence order undefined\n\n";
        }
      }
    }
  }  // end of Mult

  bool NewtonSolver::computeNewtonCorrection(
      mfem::Vector &c,
      const mfem::Vector &r,
      const mfem::Vector &u) const noexcept {
    auto profiler = this->ctx_ptr != nullptr
                        ? this->ctx_ptr->startNewProfiling(
                              "NS::computeNewtonCorrection",
                              this->ctx_ptr->isProfilingEnabled())
                        : mgis::ProfilingSection{};
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    MFEM_ASSERT(this->prec != nullptr,
                "the Solver is not set (use setLinearSolver).");
    this->prec->SetOperator(this->getJacobian(u));
    {
      auto profiler_mfem =
          this->ctx_ptr != nullptr
              ? this->ctx_ptr->startNewProfiling(
                    "MFEM::Mult(r,c)", this->ctx_ptr->isProfilingEnabled())
              : mgis::ProfilingSection{};
      this->prec->Mult(r, c);  // c = [DF(x_i)]^{-1} [F(x_i)-b]
    }
    if (this->discardLinearSolverFailure) {
      return true;
    }
    return hasConverged(*(this->prec));
  }  // end of computeNewtonCorrection

  NewtonSolver::~NewtonSolver() = default;

}  // end of namespace mfem_mgis
