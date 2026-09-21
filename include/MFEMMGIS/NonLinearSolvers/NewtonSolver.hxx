/*!
 * \file   MFEMMGIS/NonLinearSolvers/NewtonSolver.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX
#define LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX

#include <vector>
#include <optional>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverBase.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  //! \brief custom implementation of the Newton Solver
  struct MFEM_MGIS_EXPORT NewtonSolver : public NonLinearSolverBase {
#ifdef MFEM_USE_MPI
    //! \brief default constructor
    NewtonSolver(NonLinearEvolutionProblemImplementation<true> &);
#endif /* MFEM_USE_MPI */
    //! \brief default constructor
    NewtonSolver(NonLinearEvolutionProblemImplementation<false> &);
    //
    void Mult(const mfem::Vector &, mfem::Vector &) const override;
    //! \brief destructor
    ~NewtonSolver() override;

   protected:
    /*!
     * \brief compute the correction associated with the given residual
     * \param[in] c: Newton' correction
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     */
    virtual bool computeNewtonCorrection(mfem::Vector &,
                                         const mfem::Vector &,
                                         const mfem::Vector &) const noexcept;
  };  // end of struct NewtonSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX */
