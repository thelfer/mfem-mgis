/*!
 * \file   src/NonLinearSolverBase.cxx
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
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"

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

  NonLinearSolverBase::NonLinearSolverBase(
      NonLinearEvolutionProblemImplementation<true> &p)
      : AbstractNonLinearSolver(p) {
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions([&p](const mfem::Vector &u) {
      return p.integrate(
          u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR, {});
    });
  }  // end of NonLinearSolverBase

#endif /* MFEM_USE_MPI */

  NonLinearSolverBase::NonLinearSolverBase(
      NonLinearEvolutionProblemImplementation<false> &p)
      : AbstractNonLinearSolver(p) {
    checkSolverOperator(p);
    this->oper = &p;
    this->height = p.Height();
    this->width = p.Width();
    this->iterative_mode = true;
    this->addNewUnknownsEstimateActions([&p](const mfem::Vector &u) {
      return p.integrate(
          u, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR, {});
    });
  }  // end of NonLinearSolverBase

  bool NonLinearSolverBase::setSolverParameters(
      Context &ctx, const Parameters &params) noexcept {
    auto allowed_parameters = getIterativeSolverParametersList();
    allowed_parameters.push_back("DiscardLinearSolverFailure");
    if (!checkParameters(ctx, params, allowed_parameters)) {
      return false;
    }
    const auto osubparams =
        extract(ctx, params, getIterativeSolverParametersList());
    if (isInvalid(osubparams)) {
      return false;
    }
    if (!mfem_mgis::setSolverParameters(ctx, *this, *osubparams)) {
      return false;
    }
    if (contains(params, "DiscardLinearSolverFailure")) {
      const auto ob = get<bool>(ctx, params, "DiscardLinearSolverFailure");
      if (isInvalid(ob)) {
        return false;
      }
      this->discardLinearSolverFailure = *ob;
    }
    return true;
  }  // end of setSolverParameters

  void NonLinearSolverBase::SetOperator(const mfem::Operator &) {
    raise("NonLinearSolverBase::SetOperator: invalid call");
  }  // end of SetOperator

  void NonLinearSolverBase::SetPreconditioner(Solver &) {
    raise("NonLinearSolverBase::SetOperator: invalid call");
  }  // end of SetPreconditioner

  void NonLinearSolverBase::setLinearSolver(LinearSolver &s) noexcept {
    this->prec = &s;
    this->prec->iterative_mode = false;
  }  // end of setLinearSolver

  real NonLinearSolverBase::GetInitialNorm() const noexcept {
    if (this->reference_residual_norm.has_value()) {
      return *(this->reference_residual_norm);
    }
    return real{};
  }  // end of GetInitialNorm

  bool NonLinearSolverBase::setReferenceResidualNorm(Context &ctx,
                                                     const real v) noexcept {
    if (v <= 0) {
      return ctx.registerErrorMessage(
          "negative value given for the reference norm of the residual");
    }
    this->reference_residual_norm = v;
    return true;
  }  // end of setReferenceResidualNorm

  void NonLinearSolverBase::unsetReferenceResidualNorm() noexcept {
    this->reference_residual_norm.reset();
  }  // end of unsetReferenceResidualNorm

  void NonLinearSolverBase::setContext(Context &ctx) noexcept {
    this->ctx_ptr = &ctx;
  }  // end of setContext

  void NonLinearSolverBase::unsetContext() noexcept { this->ctx_ptr = nullptr; }

  std::vector<Parameter> NonLinearSolverBase::getIterationsInformation()
      const noexcept {
    return this->iterations_information;
  }  // end of getIterationsInformation

  bool NonLinearSolverBase::isLinearSolverFailureDiscarded() const noexcept {
    return this->discardLinearSolverFailure;
  }  // end of isLinearSolverFailureDiscarded

  void NonLinearSolverBase::computeResidual(mfem::Vector &r,
                                            const mfem::Vector &u) const {
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    this->oper->Mult(u, r);
  }  // end of computeResidual

  mfem::Operator &NonLinearSolverBase::getJacobian(
      const mfem::Vector &u) const {
    auto profiler =
        this->ctx_ptr != nullptr
            ? this->ctx_ptr->startNewProfiling(
                  "NS::getJacobian", this->ctx_ptr->isProfilingEnabled())
            : mgis::ProfilingSection{};
    MFEM_ASSERT(this->oper != nullptr,
                "the Operator is not set (use SetOperator).");
    return this->oper->GetGradient(u);
  }  // end of getJacobian

  void NonLinearSolverBase::addNewUnknownsEstimateActions(
      std::function<bool(const mfem::Vector &)> a) noexcept {
    auto profiler = this->ctx_ptr != nullptr
                        ? this->ctx_ptr->startNewProfiling(
                              "NS::addNewUnknownsEstimateActions",
                              this->ctx_ptr->isProfilingEnabled())
                        : mgis::ProfilingSection{};
    if (a) {
      this->nue_actions.push_back(std::move(a));
    }
  }  // end of addNewUnknownsEstimateActions

  bool NonLinearSolverBase::processNewUnknownsEstimate(
      const mfem::Vector &u) const {
    auto profiler = this->ctx_ptr != nullptr
                        ? this->ctx_ptr->startNewProfiling(
                              "NS::processNewUnknownsEstimate",
                              this->ctx_ptr->isProfilingEnabled())
                        : mgis::ProfilingSection{};
    for (const auto &a : this->nue_actions) {
      if (!a(u)) {
        return false;
      }
    }
    return true;
  }  // end of processNewUnknownsEstimate

  NonLinearSolverBase::~NonLinearSolverBase() = default;

}  // end of namespace mfem_mgis
