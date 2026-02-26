/*!
 * \file   include/MFEMMGIS/LinearSolverHandler.hxx
 * \brief
 * \author Thomas Helfer
 * \date   20/01/2026
 */

#ifndef LIB_MFEM_MGIS_LINEARSOLVERHANDLER_HXX
#define LIB_MFEM_MGIS_LINEARSOLVERHANDLER_HXX

#include "mfem/linalg/solvers.hpp"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief result from the linear solver factories
   */
  struct [[nodiscard]] LinearSolverHandler {
    std::unique_ptr<LinearSolver> linear_solver;
    std::unique_ptr<LinearSolverPreconditioner> preconditioner;
  };  // end of LinearSolverHandler

  MFEM_MGIS_EXPORT [[nodiscard]] bool isInvalid(
      const LinearSolverHandler&) noexcept;

}  // end of namespace mfem_mgis

namespace mgis::internal {

  //! \brief partial specialisation for linear solver handlers
  template <>
  struct InvalidValueTraits<mfem_mgis::LinearSolverHandler> {
    static constexpr bool isSpecialized = true;
    static mfem_mgis::LinearSolverHandler getValue() noexcept { return {}; }
  };

}  // end of namespace mgis::internal

#endif /* LIB_MFEM_MGIS_LINEARSOLVERHANDLER_HXX */
