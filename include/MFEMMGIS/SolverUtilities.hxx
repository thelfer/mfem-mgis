/*!
 * \file   MFEMMGIS/SolverUtilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   30/03/2021
 */

#ifndef LIB_MFEM_MGIS_SOLVERUTILITIES_HXX
#define LIB_MFEM_MGIS_SOLVERUTILITIES_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief return the list of parameters for an interative solver
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
   * \param[in] s: iterative solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context&, IterativeSolver&, const Parameters&) noexcept;

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in, out] ctx: execution context
   * \param[in] s: solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool setSolverParameters(
      Context&, mfem::PetscNonlinearSolver&, const Parameters&) noexcept;
#endif /* MFEM_USE_PETSC */

  /*!
   * \brief set the parameters of an iterative solver
   *
   * \param[in] throwing: dummy argument to indicate that this function may
   * throw an exception.
   * \param[in] s: iterative solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(attributes::Throwing,
                                                           IterativeSolver&,
                                                           const Parameters&);

#ifdef MFEM_USE_PETSC
  /*!
   * \brief set the parameters of a PETSc solver
   *
   * \param[in] throwing: dummy argument to indicate that this function may
   * throw an exception.
   * \param[in] s: solver
   * \param[in] params: parameters
   */
  MFEM_MGIS_EXPORT [[deprecated]] void setSolverParameters(
      attributes::Throwing, mfem::PetscNonlinearSolver&, const Parameters&);
#endif /* MFEM_USE_PETSC */

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SOLVERUTILITIES_HXX */
