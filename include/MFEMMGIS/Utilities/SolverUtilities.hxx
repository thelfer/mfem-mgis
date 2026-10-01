/*!
 * \file   MFEMMGIS/Utilities/SolverUtilities.hxx
 * \brief  This file declares functions setting the parameters of solvers and
 * checking their convergence
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX
#define LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX

#include <optional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \return the list of parameters for an iterative solver
   * \note the following parameters are available:
   * - `AbstractNonLinearEvolutionProblem::SolverVerbosityLevel`, aka
   *   `"VerbosityLevel"`,
   * - `AbstractNonLinearEvolutionProblem::SolverRelativeTolerance`, aka
   *   `"RelativeTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverAbsoluteTolerance`, aka
   *   `"AbsoluteTolerance"`,
   * - `AbstractNonLinearEvolutionProblem::SolverMaximumNumberOfIterations`, aka
   *   `"MaximumNumberOfIterations"`,
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<std::string>
  getIterativeSolverParametersList();

  /*!
   * \brief set the parameters of an iterative solver
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] s: iterative solver
   * \param[in] params: parameters
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context& ctx, IterativeSolver& s, const Parameters& params) noexcept;

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] s: solver
   * \param[in] params: parameters
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context& ctx,
      mfem::PetscNonlinearSolver& s,
      const Parameters& params) noexcept;
#endif /* MFEM_USE_PETSC */

  /*!
   * \brief set the parameters of an iterative solver
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in, out] s: iterative solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(
      attributes::Throwing throwing,
      IterativeSolver& s,
      const Parameters& params);

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in, out] s: solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(
      attributes::Throwing throwing,
      mfem::PetscNonlinearSolver& s,
      const Parameters& params);
#endif /* MFEM_USE_PETSC */

  /*!
   * \brief check if the linear solver has converged
   * \param[in] ls: linear solver
   * \return if the linear solver has converged
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool hasConverged(
      const LinearSolver& ls) noexcept;  // end of hasConverged
  /*!
   * \brief get the number of iterations of an iterative linear solver
   *
   * Unified for both LinearSolver and Hypre solvers.
   * Optional return type in case the LinearSolver is not iterative.
   *
   * \param[in] ls: linear solver
   * \return the number of iterations, empty if the solver is not iterative
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<int>
  getNumberOfIterationsAtConvergence(const LinearSolver& ls) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UTILITIES_SOLVERUTILITIES_HXX */
