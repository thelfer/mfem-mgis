/*!
 * \file   MFEMMGIS/FirstIterationConvergenceCriterion.hxx
 * \brief  This file declares the FirstIterationConvergenceCriterion class
 * \date   02/09/2024
 */

#ifndef LIB_MFEM_MGIS_FIRST_ITERATION_CONVERGENCE_CRITERION_HXX
#define LIB_MFEM_MGIS_FIRST_ITERATION_CONVERGENCE_CRITERION_HXX

#include <map>
#include "MFEMMGIS/CouplingSchemeConvergenceCriterionBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief a convergence criterion which checks that every coupling item
   * converged at the first iteration.
   */
  struct MFEM_MGIS_EXPORT FirstIterationConvergenceCriterion
      : CouplingSchemeConvergenceCriterionBase {
    //! \return a description of this criterion
    [[nodiscard]] static std::string getDescription() noexcept;
    //! \return a description of the parameters of this criterion
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \brief constructor
    FirstIterationConvergenceCriterion();
    /*!
     * \brief method called at the beginning of a time step. Does nothing.
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return true on success
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context& ctx, const TimeStep& ts) noexcept override;
    /*!
     * \brief check that the solvers of all items, including those of nested
     * coupling schemes, performed no iteration.
     * \return if the criterion is satisfied, empty on failure
     * \param[in, out] ctx: execution context
     * \param[in] o: output of all items of the coupling scheme
     */
    [[nodiscard]] std::optional<bool> check(
        Context& ctx, const ComputeNextStateOutput& o) const noexcept override;
    /*!
     * \brief update the state of the criterion. Does nothing.
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool update(Context& ctx) noexcept override;
    /*!
     * \brief revert the state of the criterion. Does nothing.
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool revert(Context& ctx) noexcept override;
    //! \brief destructor
    ~FirstIterationConvergenceCriterion() noexcept override;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FIRST_ITERATION_CONVERGENCE_CRITERION_HXX */
