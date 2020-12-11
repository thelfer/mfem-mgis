/*!
 * \file   include/MFEMMGIS/ResidualOperator.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_RESIDUALOPERATOR_HXX
#define LIB_MFEM_MGIS_RESIDUALOPERATOR_HXX

#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/operator.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/nonlinearform.hpp"

namespace mfem_mgis {

  //! \brief
  struct NonLinearEvolutionProblem;

  /*!
   * \brief class representing the residual of a non linear problem.
   * In practice, just a proxy class that delegates every computations
   * to the underlying non linear problem
   */
  struct ResidualOperator : public mfem::Operator {
    /*!
     * \brief constructor
     * \param[in] p: underlying non linear problem
     */
    ResidualOperator(NonLinearEvolutionProblem&);
    //! \brief compute the stiffness matrix
    mfem::Operator &GetGradient(const mfem::Vector &xp) const override;
    //! \brief compute the residual
    void Mult(const mfem::Vector &k, mfem::Vector &y) const override;
    //! \brief destructor
    ~ResidualOperator() override;

   private:
    //! \brief non linear problem to be solved over one time step
    NonLinearEvolutionProblem &problem;
  };  // end of struct ResidualOperator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_RESIDUALOPERATOR_HXX */
