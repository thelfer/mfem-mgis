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

  /*!
   * \brief custom implementation of the Newton Solver
   */
  struct NewtonSolver : public mfem::IterativeSolver {
#ifdef MFEM_USE_MPI
    //! \brief default constructor
    NewtonSolver(NonLinearEvolutionProblemImplementation<true> &);
#endif /* MFEM_USE_MPI */

    //! \brief default constructor
    NewtonSolver(NonLinearEvolutionProblemImplementation<false> &);

    /*!
     * \brief set the linear solver for inverting the Jacobian.
     * \param[in] s: linear solver
     */
    virtual void setLinearSolver(LinearSolver &);
    /*!
     * \brief add a new action called when a new estimate of the unknowns is
     * available.
     * \param[in] a: action
     */
    virtual void addNewUnknownsEstimateActions(
        std::function<bool(const mfem::Vector &)>);
    /*!
     * \brief compute the correction associated with the given residual
     * \param[in] c: Newton' correction
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     */
    bool computeNewtonCorrection(mfem::Vector &,
                                 const mfem::Vector &,
                                 const mfem::Vector &) const;
    /*!
     * \brief compute the residual
     * \param[in] r: residual
     * \param[in] u: current estimate of the unknowns
     */
    void computeResidual(mfem::Vector &, const mfem::Vector &) const;
    //
    [[noreturn]] void SetPreconditioner(Solver &) override;
    [[noreturn]] void SetOperator(const mfem::Operator &) override;
    void Mult(const mfem::Vector &, mfem::Vector &) const override;

    //! \brief destructor
    ~NewtonSolver() override;

   protected:
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
