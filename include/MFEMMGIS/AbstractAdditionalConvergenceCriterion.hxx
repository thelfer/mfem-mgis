/*!
 * \file   MFEMMGIS/AbstractAdditionalConvergenceCriterion.hxx
 * \brief
 * \author Nicolas Ségala
 * \date   6/08/2026
 */

#ifndef LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX
#define LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
namespace mfem_mgis::nonlinear_solver {
    /*!
    * \brief an abstract class for additional convergence criteria 
    *
    *
    */
    struct MFEM_MGIS_EXPORT AbstractAdditionalConvergenceCriterion{
        //! \brief custom struct for the actions taken when checking if the non-linear solver has converged
        struct CheckArguments {
            const real residual_norm; //!< current norm of the residual
            const real reference_residual_norm; //!< norm of the residual obtained in the prediction or first iteration
            const int iter; //!< current number of iterations of the non-linear solver
            const int max_iter; //!< maximum number of iterations allowed for the non-linear solver
            const bool converged; //!< current convergence status of the non-linear solver
            const mfem::Vector &u; //!< current estimate of the unknowns
        };  // end of struct CheckArguments
        
        /*!
         * \brief helper function
         *
         * Called in `NonLinearEvolutionProblemImplementationBase::setup`, via `NewtonSolver::processAdditionalConvergenceCriterionHelper`. 
         * Can be used to ensure some parameters are set correctly for the prediction. 
         */
        virtual void helper() noexcept = 0;
        /*!
         * \brief reset function
         *
         * Called in `NonLinearEvolutionProblemImplementationBase::solve`, via `NewtonSolver::processAdditionalConvergenceCriterionReset`, 
         * after the prediction is computed. 
         * Can be used to set parameters to another value for the non-linear solver iterations. 
         */
        virtual void reset() noexcept = 0;
        /*!
         * \brief check function
         *
         * Checks the convergence of the non-linear solver according to some custom criterion.
         * Called in `NewtonSolver::Mult` via `NewtonSolver::processAdditionalConvergenceCriterionCheck`. 
         * Can be used to manipulate the value of parameters until a target is achieved.
         */
        [[nodiscard]] virtual std::optional<bool> check(Context&, const CheckArguments&) noexcept = 0;
        virtual ~AbstractAdditionalConvergenceCriterion() = default;
    };  // end of struct AbstractAdditionalConvergenceCriterion

}  // end of namespace mfem_mgis::nonlinear_solver

#endif /* LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX */
