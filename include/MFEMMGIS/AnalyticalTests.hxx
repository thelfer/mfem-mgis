/*!
 * \file   include/MFEMMGIS/AnalyticalTests.hxx
 * \brief
 * \author Thomas Helfer
 * \date   25/03/2021
 */

#ifndef LIB_MFEM_MGIS_ANALYTICALTESTS_HXX
#define LIB_MFEM_MGIS_ANALYTICALTESTS_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"

namespace mfem_mgis {

  // forward declaration
  struct NonLinearEvolutionProblem;

  /*!
   * \brief compute the L2 error against an analytical solution
   * \return the L2 norm of the difference between the unknowns at the end of
   * the time step and an analytical solution.
   * \param[in] p: considered problem
   * \param[in] f: reference function to compare with. It sets its first
   * argument to the solution at the point given by its second argument.
   */
  MFEM_MGIS_EXPORT real computeL2ErrorAgainstAnalyticalSolution(
      NonLinearEvolutionProblem& p,
      std::function<void(mfem::Vector&, const mfem::Vector&)> f) noexcept;

  /*!
   * \brief compare the results to an analytical solution with a specified
   * threshold using the L2 norm.
   * \return whether the L2 error is below the threshold, empty on failure
   * \param[in, out] ctx: execution context
   * \param[in] p: considered problem
   * \param[in] f: reference function to compare with
   * \param[in] params: parameters. `CriterionThreshold` is required,
   * `VerbosityLevel` is optional.
   */
  MFEM_MGIS_EXPORT std::optional<bool> compareToAnalyticalSolution(
      Context& ctx,
      NonLinearEvolutionProblem& p,
      std::function<void(mfem::Vector&, const mfem::Vector&)> f,
      const Parameters& params) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ANALYTICALTESTS_HXX */
