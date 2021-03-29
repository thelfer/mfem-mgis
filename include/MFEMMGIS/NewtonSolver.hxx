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

namespace mfem_mgis {

  /// Newton's method for solving F(x)=b for a given operator F.
  /** The method GetGradient() must be implemented for the operator F.
      The preconditioner is used (in non-iterative mode) to evaluate
      the action of the inverse gradient of the operator. */
  struct MFEM_MGIS_EXPORT NewtonSolver : public mfem::IterativeSolver {
    //! \brief default constructor
    NewtonSolver();

#ifdef MFEM_USE_MPI
    /*!
     * \brief default constructor
     * \param[in] c: communicator
     */
    NewtonSolver(MPI_Comm);
#endif

    //
    void SetOperator(const mfem::Operator &op) override;
    void Mult(const mfem::Vector &b, mfem::Vector &x) const override;
    //
    /*!
     * \brief set the linear solver for inverting the Jacobian.
     * \param[in] s: linear solver
     * \note this method is equivalent to calling SetPreconditioner().
     */
    virtual void setLinearSolver(LinearSolver &);
    /*!
     * \brief add a new action called when a new estimate of the unknowns is
     * available.
     * \param[in] a: action
     */
    virtual void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)>);
    //! \brief destructor
    ~NewtonSolver() override;

   protected:
    /** @brief This method can be overloaded in derived classes to implement
       line search algorithms. */
    /** The base class implementation (NewtonSolver) simply returns 1. A return
        value of 0 indicates a failure, interrupting the Newton iteration. */
    virtual double ComputeScalingFactor(const mfem::Vector &,
                                        const mfem::Vector &) const {
      return 1.0;
    }

    /*!
     * \brief method called when a new estimate of the unknowns is available.
     * \param[in] u: new unknown estimate
     */
    virtual bool processNewUnknownsEstimate(const mfem::Vector &) const;
    /*!
     * \brief actions performed when a new estimate of the unknowns are
     * available
     */
    std::vector<std::function<bool(const mfem::Vector &)>> nue_actions;
  };  // end of struct NewtonSolver

}  // end of namespace mfem_mgis

#endif /* LIB_NEWTONSOLVER_HXX */
