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

  /*!
   * \brief a coupling scheme calling all coupling items until all convergence
   * criteria are satisfied
   */
  struct MFEM_MGIS_EXPORT IterativeCouplingScheme : CouplingSchemeBase {
    //! \return a description of this scheme
    static std::string getDescription() noexcept;
    //! \return a description of the parameters of this scheme
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    IterativeCouplingScheme(Context &ctx, const MeshDiscretization &m);
    /*!
     * \brief set the maximum number of iterations
     * \param[in, out] ctx: execution context
     * \param[in] n: maximum number of iterations
     * \return true on success
     */
    [[nodiscard]] bool setMaximumNumberOfIterations(Context &ctx,
                                                    const size_type n) noexcept;
    //! \return the name of the scheme, `IterativeCouplingScheme` by default
    [[nodiscard]] std::string getName() const noexcept override;
    /*!
     * \return a description of the coupling scheme
     * \param[in, out] ctx: execution context
     * \param[in] b: boolean being the default value for information requests
     * \param[in] parameters: information requests. Supported requests are
     * `ShortDescription`, `NumericalParameters` and `CouplingItems`.
     */
    [[nodiscard]] std::optional<std::string> describe(
        Context &ctx,
        const bool b,
        const Parameters &parameters) const noexcept override;
    //     [[nodiscard]] bool addConvergenceCriterion(
    //         Context &, std::string_view, const Parameters &) noexcept
    //         override;
    /*!
     * \brief add a new convergence criterion
     * \param[in, out] ctx: execution context
     * \param[in] c: convergence criterion
     * \return true on success
     */
    [[nodiscard]] bool addConvergenceCriterion(
        Context &ctx,
        std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion> c) noexcept
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
    /*!
     * \brief call `performInitializationTaksAtTheBeginningOfTheTimeStep` on
     * all coupling items and all convergence criteria
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return true on success
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief call all coupling items until all convergence criteria are
     * satisfied or the maximum number of iterations is reached
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return the exit status and, on success, the outputs of the scheme
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief update all coupling items and all convergence criteria
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    /*!
     * \brief revert all coupling items and all convergence criteria
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~IterativeCouplingScheme() noexcept override;

   private:
    //! \brief list of convergence criteria
    std::vector<std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>>
        convergence_criteria;
    //! \brief maximum number of iterations
    size_type maximum_number_of_iterations = 1;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX */
