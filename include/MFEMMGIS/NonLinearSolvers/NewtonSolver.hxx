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
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NewtonSolver(NonLinearEvolutionProblemImplementation<true>& p);
#endif /* MFEM_USE_MPI */
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     */
    NewtonSolver(NonLinearEvolutionProblemImplementation<false>& p);
    /*!
     * \brief solve the non linear problem
     * \param[in] b: right hand side, unused
     * \param[in, out] x: initial guess, then solution
     */
    void Mult(const mfem::Vector& b, mfem::Vector& x) const override;
    //! \brief destructor
    ~NewtonSolver() override;

   protected:
    /*!
     * \brief compute the correction associated with the given residual
     * \param[out] c: opposite of the Newton correction
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     * \return true on success
     */
    virtual bool computeNewtonCorrection(mfem::Vector& c,
                                         const mfem::Vector& r,
                                         const mfem::Vector& u) const noexcept;
  };  // end of struct NewtonSolver

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NEWTONSOLVER_HXX */
