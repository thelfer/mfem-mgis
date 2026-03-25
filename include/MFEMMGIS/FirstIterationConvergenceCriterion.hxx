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
   * \brief a convergence criterion which checks that every coupling items
   * converged at the first iteration.
   */
  struct FirstIterationConvergenceCriterion
      : CouplingSchemeConvergenceCriterionBase {
    //! \return a description of this criterion
    [[nodiscard]] static std::string getDescription() noexcept;
    //! \return a description of the parameters of this criterion
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    //! \brief constructor
    FirstIterationConvergenceCriterion();
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept override;
    [[nodiscard]] std::optional<bool> check(
        Context &, const ComputeNextStateOutput &) const noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
    //! \brief destructor
    ~FirstIterationConvergenceCriterion() noexcept override;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FIRST_ITERATION_CONVERGENCE_CRITERION_HXX */
