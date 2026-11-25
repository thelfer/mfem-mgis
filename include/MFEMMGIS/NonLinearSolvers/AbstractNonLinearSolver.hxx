/*!
 * \file   MFEMMGIS/NonLinearSolvers/AbstractNonLinearSolver.hxx
 * \brief  This file declares the `AbstractNonLinearSolver` class
 * \author Thomas Helfer
 * \date   20/09/2026
 */

#ifndef LIB_MFEMMGIS_NONLINEARSOLVERS_ABSTRACTNONLINEARSOLVER_HXX
#define LIB_MFEMMGIS_NONLINEARSOLVERS_ABSTRACTNONLINEARSOLVER_HXX

#include <memory>
#include <vector>
#include <optional>
#include "mfem/linalg/solvers.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declarations
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief base class for nonlinear solvers
   *
   * \note This class inherits from `mfem::IterativeSolver`. It is thus not
   * purely virtual. Child classes must implement `mfem::IterativeSolver::Mult`
   */
  struct MFEM_MGIS_EXPORT AbstractNonLinearSolver
      : public mfem::IterativeSolver {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    AbstractNonLinearSolver(NonLinearEvolutionProblemImplementation<true>& p);
#endif /* MFEM_USE_MPI */
    //! \brief constructor
    AbstractNonLinearSolver(NonLinearEvolutionProblemImplementation<false>&);
    /*!
     * \brief set the solver parameters
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters
     * \return true on success
     */
    [[nodiscard]] virtual bool setSolverParameters(
        Context& ctx, const Parameters& params) noexcept = 0;
    //! \return if the failure of the linear solver is discarded
    [[nodiscard]] virtual bool isLinearSolverFailureDiscarded()
        const noexcept = 0;
    /*!
     * \brief set the linear solver for inverting the Jacobian.
     * \param[in] s: linear solver
     */
    virtual void setLinearSolver(LinearSolver& s) noexcept = 0;
    /*!
     * \brief add a new action called when a new estimate of the unknowns is
     * available.
     * \param[in] a: action
     */
    virtual void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)> a) noexcept = 0;
    /*!
     * \brief register an additional convergence criterion
     * \param[in] a: additional convergence criterion
     */
    virtual void addAdditionalConvergenceCriterion(
        std::shared_ptr<AbstractAdditionalConvergenceCriterion> a) = 0;
    /*!
     * \brief evaluate all registered additional convergence criteria
     * \param[in, out] ctx: execution context
     * \param[in] s: arguments passed to
     * `AbstractAdditionalConvergenceCriterion::check`
     * \return the combined convergence status, or `std::nullopt` if a criterion
     * reported an error
     */
    virtual std::optional<bool> processAdditionalConvergenceCriterionCheck(
        Context &ctx,
        const AbstractAdditionalConvergenceCriterion::CheckArguments &s)
        const = 0;
    /*!
     * \brief method called after the prediction step and before the nonlinear
     * resolution, see `NonLinearEvolutionProblemImplementationBase::solve`
     */
    virtual void processAdditionalConvergenceCriterionReset() = 0;
    /*!
     * \brief method called during problem setup before the prediction step,
     * see `NonLinearEvolutionProblemImplementationBase::setup`
     */
    virtual void processAdditionalConvergenceCriterionHelper() = 0;
    /*!
     * \brief set the reference value for the norm of the residual.
     * \param[in, out] ctx: execution context
     * \param[in] v: value of the reference residual
     * \return true on success
     */
    [[nodiscard]] virtual bool setReferenceResidualNorm(
        Context& ctx, const real v) noexcept = 0;
    //! \brief unset the reference residual norm
    virtual void unsetReferenceResidualNorm() noexcept = 0;
    //! \return the initial norm of the residual
    [[nodiscard]] virtual real GetInitialNorm() const noexcept = 0;
    /*!
     * \brief set the current execution context
     *
     * \param[in, out] ctx: execution context
     *
     * \note this execution context is used inside the `Mult` method whose
     * prototype is imposed by `MFEM`'s API.
     */
    virtual void setContext(Context& ctx) noexcept = 0;
    //! \brief unset the execution context
    virtual void unsetContext() noexcept = 0;
    //! \return the information collected during the iterations
    virtual std::vector<Parameter> getIterationsInformation()
        const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearSolver() noexcept;
  };  // end of AbstractNonLinearSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARSOLVERS_ABSTRACTNONLINEARSOLVER_HXX */
