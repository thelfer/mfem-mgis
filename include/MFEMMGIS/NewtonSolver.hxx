/*!
 * \file   include/MFEMMGIS/NewtonSolver.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2021
 */

#ifndef LIB_MFEM_MGIS_NEWTONSOLVER_HXX
#define LIB_MFEM_MGIS_NEWTONSOLVER_HXX

#include <vector>
#include <functional>
#include "mfem/linalg/solvers.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  /// Newton's method for solving F(x)=b for a given operator F.
  /** The method GetGradient() must be implemented for the operator F.
      The preconditioner is used (in non-iterative mode) to evaluate
      the action of the inverse gradient of the operator. */
  template <bool parallel>
  struct NewtonSolver : public mfem::IterativeSolver {
    //! \brief default constructor
    NewtonSolver(NonLinearEvolutionProblemImplementation<parallel> &);
    /*!
     * \brief set the linear solver for inverting the Jacobian.
     * \param[in] s: linear solver
     */
    virtual void setLinearSolver(std::unique_ptr<LinearSolver>);
    //! \brief solve the problem.
    virtual void solve() const;
    /*!
     * \brief add a new action called when a new estimate of the unknowns is
     * available.
     * \param[in] a: action
     */
    virtual void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)>);
    //
    [[noreturn]] void SetPreconditioner(Solver &) override;
    [[noreturn]] void SetOperator(const mfem::Operator &) override;
    [[noreturn]] void Mult(const mfem::Vector &, mfem::Vector &) const override;
    //! \brief destructor
    ~NewtonSolver() override;

   protected:

    /*!
     * \brief method called when a new estimate of the unknowns is available.
     * \param[in] u: new unknown estimate
     */
    virtual bool processNewUnknownsEstimate(const mfem::Vector &) const;
    //! \brief underlying problem
    NonLinearEvolutionProblemImplementation<parallel> &problem;
    /*!
     * \brief actions performed when a new estimate of the unknowns are
     * available
     */
    std::vector<std::function<bool(const mfem::Vector &)>> nue_actions;
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    //!
    bool prediction = true;
  };  // end of struct NewtonSolver


#ifdef MFEM_USE_MPI

  template <>
  NewtonSolver<true>::NewtonSolver(
      NonLinearEvolutionProblemImplementation<true> &);

#endif MFEM_USE_MPI

  template <>
  NewtonSolver<false>::NewtonSolver(
      NonLinearEvolutionProblemImplementation<false> &);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/NewtonSolver.ixx"

#endif /* LIB_NEWTONSOLVER_HXX */
