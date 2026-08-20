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
    
    

    struct MFEM_MGIS_EXPORT AbstractAdditionalConvergenceCriterion{
        //! \brief custom struct for the actions after the non-linear solver has converged
        /* \param[in] residual_norm TODO
        * \param[in] reference_residual_norm TODO
        * \param[in] iterations TODO
        * \param[in] max_iterations TODO
        * \param[in] convergence TODO
        * \param[in] u: current estimate of the unknowns
        */
        struct CheckArguments {
            const real residual_norm;
            const real reference_residual_norm;
            const int iter;
            const int max_iter;
            const bool converged;
            const mfem::Vector &u;
        };  // end of struct CheckArguments
        


        virtual void reset() noexcept = 0;
        [[nodiscard]] virtual std::optional<bool> check(Context&, const CheckArguments&) noexcept = 0;
        virtual ~AbstractAdditionalConvergenceCriterion() = default;
    };  // end of struct AbstractAdditionalConvergenceCriterion

}  // end of namespace mfem_mgis::nonlinear_solver

#endif /* LIB_MFEM_MGIS_ABSTRACTADDITIONALCONVERGENCECRITERION_HXX */
