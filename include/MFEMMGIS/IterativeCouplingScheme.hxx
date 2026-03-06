/*!
 * \file   MFEMMGIS/IterativeCouplingScheme.hxx
 * \brief  This file declares the `IterativeCouplingScheme` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_ITERATIVE_COUPLING_SCHEME_HXX
#define LIB_MFEM_MGIS_ITERATIVE_COUPLING_SCHEME_HXX

#include <map>
#include <string>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/CouplingSchemeBase.hxx"

namespace mfem_mgis {

  /*!´
   * \brief the simpliest coupling scheme: all declared models are called once
   */
  struct MFEM_MGIS_EXPORT IterativeCouplingScheme : CouplingSchemeBase {
    //! \return a description of this scheme
    static std::string getDescription() noexcept;
    //! \return a description of the parameters of this scheme
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    IterativeCouplingScheme(const MeshDiscretization &);
    /*!
     * \brief set the maximum number of iterations
     * \param[in] ctx: execution context
     * \param[in] n: maximum number of iterations
     */
    [[nodiscard]] bool setMaximumNumberOfIterations(Context &,
                                                    const size_type) noexcept;
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] std::optional<std::string> describe(
        Context &, const bool, const Parameters &) const noexcept override;
    //     [[nodiscard]] bool addConvergenceCriterion(
    //         Context &, std::string_view, const Parameters &) noexcept
    //         override;
    [[nodiscard]] bool addConvergenceCriterion(
        Context &,
        std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>) noexcept
        override;
    //     [[nodiscard]] bool declareDependencies(
    //         Context &, DependenciesManager &) const noexcept override;
    //     [[nodiscard]] bool initializeBeforeResourcesAllocation(
    //         Context &,
    //         ValueEvaluatorsFactory &,
    //         NodalEvaluatorsFactory &,
    //         IPEvaluatorsFactory &) noexcept override;
    //     [[nodiscard]] bool initializeAfterResourcesAllocation(Context &)
    //     noexcept override;
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept override;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
    //! \brief destructor
    ~IterativeCouplingScheme() noexcept override;

   private:
    //! \brief list of convergence criteria
    std::vector<std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>>
        convergence_criteria;
    //! \brief number of iterations
    size_type maximum_number_of_iterations = 1;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX */
