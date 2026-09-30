/*!
 * \file   include/MFEMMGIS/LinearSolverHandler.hxx
 * \brief  This file declares the `LinearSolverHandler` class
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
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    //! \brief preconditioner
    std::unique_ptr<LinearSolverPreconditioner> preconditioner;
  };  // end of LinearSolverHandler

  /*!
   * \brief check if the given linear solver handler is invalid
   * \param[in] s: linear solver handler
   * \return if the given linear solver handler is invalid
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isInvalid(
      const LinearSolverHandler& s) noexcept;

}  // end of namespace mfem_mgis

namespace mgis::internal {

  //! \brief specialisation for linear solver handlers
  template <>
  struct InvalidValueTraits<mfem_mgis::LinearSolverHandler> {
    //! \brief tag indicating that this class is properly specialized
    static constexpr bool isSpecialized = true;
    //! \return an invalid linear solver handler
    static mfem_mgis::LinearSolverHandler getValue() noexcept { return {}; }
  };

}  // end of namespace mgis::internal

#endif /* LIB_MFEM_MGIS_LINEARSOLVERHANDLER_HXX */
