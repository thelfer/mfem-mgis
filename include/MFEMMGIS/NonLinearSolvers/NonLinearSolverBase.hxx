/*!
 * \file   MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx
 * \brief  This file declares the `NonLinearSolverBase` class
 * \author Thomas Helfer
 * \date   20/09/2026
 */

#ifndef LIB_MFEMMGIS_NONLINEARSOLVERS_NONLINEARSOLVERBASE_HXX
#define LIB_MFEMMGIS_NONLINEARSOLVERS_NONLINEARSOLVERBASE_HXX

#include "mfem/linalg/solvers.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/AbstractNonLinearSolver.hxx"

namespace mfem_mgis {

  //! \brief base class for nonlinear solvers
  struct MFEM_MGIS_EXPORT NonLinearSolverBase : public AbstractNonLinearSolver {
#ifdef MFEM_USE_MPI
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NonLinearSolverBase(NonLinearEvolutionProblemImplementation<true> &p);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NonLinearSolverBase(NonLinearEvolutionProblemImplementation<false> &p);
    //
    void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)> a) noexcept override final;
    //
    [[nodiscard]] bool setSolverParameters(
        Context &ctx, const Parameters &params) noexcept override;
    [[nodiscard]] bool isLinearSolverFailureDiscarded() const noexcept override;
    void setLinearSolver(LinearSolver &s) noexcept override;
    /*!
     * \brief set the reference value for the norm of the residual.
     * \param[in, out] ctx: execution context
     * \param[in] v: value of the reference residual, must be positive
     * \return true on success
     */
    [[nodiscard]] bool setReferenceResidualNorm(Context &ctx,
                                                const real v) noexcept override;
    void unsetReferenceResidualNorm() noexcept override;
    //! \return the reference residual norm if set, zero otherwise
    [[nodiscard]] real GetInitialNorm() const noexcept override;
    void setContext(Context &ctx) noexcept override;
    void unsetContext() noexcept override;
    std::vector<Parameter> getIterationsInformation() const noexcept override;
    //! \brief destructor
    ~NonLinearSolverBase() noexcept;

   protected:
    /*!
     * \brief not supported, throws an exception
     * \param[in] s: preconditioner, unused
     */
    [[noreturn]] void SetPreconditioner(Solver &s) override;
    /*!
     * \brief not supported, throws an exception
     * \param[in] op: operator, unused
     */
    [[noreturn]] void SetOperator(const mfem::Operator &op) override;
    /*!
     * \brief compute the residual
     * \param[out] r: residual
     * \param[in] u: current estimate of the unknowns
     */
    virtual void computeResidual(mfem::Vector &r, const mfem::Vector &u) const;
    /*!
     * \brief return the jacobian of the system
     * \param[in] u: current estimate of the unknowns
     * \return the jacobian of the system
     */
    virtual mfem::Operator &getJacobian(const mfem::Vector &u) const;
    /*!
     * \brief method called when a new estimate of the unknowns is available.
     * \param[in] u: new unknown estimate
     * \return true on success
     */
    virtual bool processNewUnknownsEstimate(const mfem::Vector &u) const;
    /*!
     * \brief actions performed when a new estimate of the unknowns is
     * available
     */
    std::vector<std::function<bool(const mfem::Vector &)>> nue_actions;
    /*!
     * \brief information collected during the iterations
     * This vector is meant to be cleared at the beginning of the Mult method
     */
    mutable std::vector<Parameter> iterations_information;
    /*!
     * \brief data containing the reference value for the norm of the residual.
     *
     * \note This value can be set before calling `Mult` when a prediction of
     * the solution is made. If this value is not set, it is set to the value of
     * the residual at the first iteration.
     */
    mutable std::optional<real> reference_residual_norm;
    //! \brief pointer to an execution context
    Context *ctx_ptr = nullptr;
    /*!
     * \brief boolean stating if failure of the linear solver must not be
     * checked
     */
    bool discardLinearSolverFailure = false;
  };  // end of NonLinearSolverBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEARSOLVERS_NONLINEARSOLVERBASE_HXX */
