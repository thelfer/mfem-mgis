/*!
 * \file   MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx
 * \brief  This file declares the `AbstractAdditionalConvergenceCriterion` class
 * \author Nicolas Ségala
 * \date   6/08/2026
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX
#define LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {
  /*!
   * \brief an abstract class for additional convergence criteria
   */
  struct MFEM_MGIS_EXPORT AbstractAdditionalConvergenceCriterion {
    //! \brief arguments passed when checking nonlinear solver convergence
    struct CheckArguments {
      const real residual_norm;            //!< current norm of the residual
      const real reference_residual_norm;  //!< norm of the residual obtained
                                           //!< in the prediction or first
                                           //!< iteration
      const int iter;      //!< current iteration index of the nonlinear solver
      const int max_iter;  //!< maximum number of iterations allowed for the
                           //!< nonlinear solver
      const bool converged;  //!< current convergence status of the nonlinear
                             //!< solver
      const mfem::Vector &u;  //!< current estimate of the unknowns
    };                        // end of struct CheckArguments

    /*!
     * \brief helper function
     *
     * Called in `NonLinearEvolutionProblemImplementationBase::setup`, via
     * `AbstractNonLinearSolver::processAdditionalConvergenceCriterionHelper`.
     * Can be used to ensure parameters are set appropriately before the
     * prediction step.
     */
    virtual void helper() noexcept = 0;
    /*!
     * \brief reset function
     *
     * Called in `NonLinearEvolutionProblemImplementationBase::solve`, via
     * `AbstractNonLinearSolver::processAdditionalConvergenceCriterionReset`,
     * after the prediction is computed and before the nonlinear solver runs.
     * Can be used to initialize or arm parameters for the nonlinear solver
     * iterations.
     */
    virtual void reset() noexcept = 0;
    /*!
     * \brief check function
     *
     * Checks the convergence of the nonlinear solver according to a custom
     * criterion. Called in `NewtonSolver::Mult` via
     * `AbstractNonLinearSolver::processAdditionalConvergenceCriterionCheck`.
     * Can also be used to update solver, material, or behaviour parameters
     * during iterations until a target is achieved.
     *
     * \param[in, out] ctx: execution context
     * \param[in] args: convergence state and iteration arguments
     * \return `true` if the criterion is satisfied, `false` to prevent
     * convergence and continue iterating, or `std::nullopt` (`InvalidResult`)
     * to report an error and abort the nonlinear resolution.
     */
    [[nodiscard]] virtual std::optional<bool> check(
        Context &ctx, const CheckArguments &args) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractAdditionalConvergenceCriterion() noexcept = default;
  };  // end of struct AbstractAdditionalConvergenceCriterion

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX */
