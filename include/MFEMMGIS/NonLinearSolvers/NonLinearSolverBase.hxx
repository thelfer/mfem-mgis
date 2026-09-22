/*!
 * \file   MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx
 * \brief  This file declars the `NonLinearSolverBase` class
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
    //! \brief default constructor
    NonLinearSolverBase(NonLinearEvolutionProblemImplementation<true> &);
#endif /* MFEM_USE_MPI */
    //! \brief default constructor
    NonLinearSolverBase(NonLinearEvolutionProblemImplementation<false> &);
    //
    void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)>) noexcept override final;
    //
    [[nodiscard]] bool setSolverParameters(
        Context &, const Parameters &) noexcept override;
    [[nodiscard]] bool isLinearSolverFailureDiscarded() const noexcept override;
    void setLinearSolver(LinearSolver &) noexcept override;
    [[nodiscard]] bool setReferenceResidualNorm(Context &,
                                                const real) noexcept override;
    void unsetReferenceResidualNorm() noexcept override;
    [[nodiscard]] real GetInitialNorm() const noexcept override;
    void setContext(Context &) noexcept override;
    void unsetContext() noexcept override;
    std::vector<Parameter> getIterationsInformation() const noexcept;
    //! \brief destructor
    ~NonLinearSolverBase() noexcept;

   protected:
    [[noreturn]] void SetPreconditioner(Solver &) override;
    [[noreturn]] void SetOperator(const mfem::Operator &) override;
    /*!
     * \brief compute the residual
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     */
    virtual void computeResidual(mfem::Vector &, const mfem::Vector &) const;
    //! \return the jacobian of the system
    virtual mfem::Operator &getJacobian(const mfem::Vector &) const;
    /*!
     * \brief method called when a new estimate of the unknowns is available.
     * \param[in] u: new unknown estimate
     */
    virtual bool processNewUnknownsEstimate(const mfem::Vector &) const;
    /*!
     * \brief actions performed when a new estimate of the unknowns are
     * available
     */
    std::vector<std::function<bool(const mfem::Vector &)>> nue_actions;
    /*!
     * \return the information collected during the iterations
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
